// ES-EKF Hybrid Optimized for ESP32
#include <Wire.h>
#include <SparkFun_u-blox_GNSS_v3.h>
#include <esp_now.h>
#include <WiFi.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <REG.h>
#include <wit_c_sdk.h>
#include <freertos/semphr.h>

// --- CẤU HÌNH HYBRID ---
typedef float esekf_float_t;   // Float cho tính toán ESKF internal
typedef double esekf_double_t; // Double cho GNSS và chuyển đổi tọa độ

const int NOMINAL_STATE_DIM = 5;
const int ERROR_STATE_DIM = 8;

const esekf_float_t V_MIN_HEADING = 1.0f;
const esekf_float_t M_PI_F = 3.14159265358979323846f;
const esekf_double_t M_PI_D = 3.14159265358979323846;

const uint8_t SYNC_BYTE = 0xAA; // Sync header

struct ESKF_STATE_DATA
{
  esekf_float_t state[NOMINAL_STATE_DIM]; // Float cho output
  esekf_float_t P_diag[ERROR_STATE_DIM];
  esekf_float_t raw_data[9];
  esekf_float_t gnss_accuracy[3];
  uint32_t timestamp;
} __attribute__((packed));

struct TUNING_CONFIG
{
  esekf_float_t Q_noise_params[3];
  esekf_float_t R_noise_params[3];
} __attribute__((packed));

TUNING_CONFIG currentConfig;

#define RXD2 16
#define TXD2 17
#define SENSOR_SLAVE_ID 0x50

// --- BIẾN TOÀN CỤC HYBRID ---
SFE_UBLOX_GNSS myGNSS;
SemaphoreHandle_t stateMutex;

// ESKF internal states - SỬ DỤNG FLOAT
static esekf_float_t nominal_state[NOMINAL_STATE_DIM] = {0};
static esekf_float_t error_state[ERROR_STATE_DIM] = {0};
static esekf_float_t P[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
static esekf_float_t Q[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
esekf_float_t R_yaw;
esekf_float_t C0_yaw = 0.0f;

// GNSS và tọa độ - SỬ DỤNG DOUBLE
struct PointLLH
{
  esekf_double_t latitude, longitude, height;
};
struct PointECEF
{
  esekf_double_t x, y, z;
};
PointLLH initialLLH;
PointECEF initialECEF;
PointECEF enu_origin_ecef;
esekf_double_t ecef_to_enu_matrix[3][3];

const esekf_double_t WGS84_A = 6378137.0;
const esekf_double_t WGS84_E2 = 0.00669437999014;

// --- LỚP TỐI ƯU HÓA TÍNH TOÁN ---
class HybridMath
{
private:
  static const int TABLE_SIZE = 360;
  esekf_float_t sin_table[TABLE_SIZE];
  esekf_float_t cos_table[TABLE_SIZE];

public:
  HybridMath()
  {
    for (int i = 0; i < TABLE_SIZE; i++)
    {
      esekf_float_t angle_rad = i * M_PI_F / 180.0f;
      sin_table[i] = sinf(angle_rad);
      cos_table[i] = cosf(angle_rad);
    }
  }

  esekf_float_t fast_sin(esekf_float_t angle_rad)
  {
    int idx = ((int)(angle_rad * 180.0f / M_PI_F) % TABLE_SIZE + TABLE_SIZE) % TABLE_SIZE;
    return sin_table[idx];
  }

  esekf_float_t fast_cos(esekf_float_t angle_rad)
  {
    int idx = ((int)(angle_rad * 180.0f / M_PI_F) % TABLE_SIZE + TABLE_SIZE) % TABLE_SIZE;
    return cos_table[idx];
  }

  // Double precision functions for critical calculations
  esekf_double_t accurate_sin(esekf_double_t angle_rad)
  {
    return sin(angle_rad);
  }

  esekf_double_t accurate_cos(esekf_double_t angle_rad)
  {
    return cos(angle_rad);
  }
};

HybridMath hybrid_math;

// --- HÀM MA TRẬN TỐI ƯU (FLOAT) ---
void matrixMultiply5x5_5x5_float(const esekf_float_t *A, const esekf_float_t *B, esekf_float_t *C)
{
  for (int i = 0; i < 5; i++)
  {
    for (int j = 0; j < 5; j++)
    {
      C[i * 5 + j] = A[i * 5 + 0] * B[0 * 5 + j] + A[i * 5 + 1] * B[1 * 5 + j] + A[i * 5 + 2] * B[2 * 5 + j] +
                     A[i * 5 + 3] * B[3 * 5 + j] + A[i * 5 + 4] * B[4 * 5 + j];
    }
  }
}

bool inverse_5x5_optimized_float(const esekf_float_t *A, esekf_float_t *A_inv)
{
  esekf_float_t LU[25];
  memcpy(LU, A, sizeof(LU));

  int pivot[5];

  for (int i = 0; i < 5; i++)
  {
    int max_row = i;
    for (int j = i + 1; j < 5; j++)
    {
      if (fabsf(LU[j * 5 + i]) > fabsf(LU[max_row * 5 + i]))
      {
        max_row = j;
      }
    }

    if (fabsf(LU[max_row * 5 + i]) < 1e-6f)
      return false;

    if (max_row != i)
    {
      for (int j = 0; j < 5; j++)
      {
        esekf_float_t temp = LU[i * 5 + j];
        LU[i * 5 + j] = LU[max_row * 5 + j];
        LU[max_row * 5 + j] = temp;
      }
    }

    pivot[i] = max_row;

    for (int j = i + 1; j < 5; j++)
    {
      LU[j * 5 + i] /= LU[i * 5 + i];
      for (int k = i + 1; k < 5; k++)
      {
        LU[j * 5 + k] -= LU[j * 5 + i] * LU[i * 5 + k];
      }
    }
  }

  for (int col = 0; col < 5; col++)
  {
    esekf_float_t x[5] = {0};

    for (int i = 0; i < 5; i++)
    {
      x[i] = (col == pivot[i]) ? 1.0f : 0.0f;
      for (int j = 0; j < i; j++)
      {
        x[i] -= LU[i * 5 + j] * x[j];
      }
    }

    for (int i = 4; i >= 0; i--)
    {
      for (int j = i + 1; j < 5; j++)
      {
        x[i] -= LU[i * 5 + j] * x[j];
      }
      x[i] /= LU[i * 5 + i];
    }

    for (int i = 0; i < 5; i++)
    {
      A_inv[i * 5 + col] = x[i];
    }
  }

  return true;
}

esekf_float_t normalizeAngle_float(esekf_float_t angle)
{
  return atan2f(sinf(angle), cosf(angle));
}

// --- HÀM GNSS VÀ TỌA ĐỘ (DOUBLE PRECISION) ---
void GetRTK(PointLLH &_CurrentPos)
{
  int32_t latitude = myGNSS.getHighResLatitude();
  int8_t latitudeHp = myGNSS.getHighResLatitudeHp();
  _CurrentPos.latitude = ((esekf_double_t)latitude / 10000000.0) + ((esekf_double_t)latitudeHp / 1000000000.0);

  int32_t longitude = myGNSS.getHighResLongitude();
  int8_t longitudeHp = myGNSS.getHighResLongitudeHp();
  _CurrentPos.longitude = ((esekf_double_t)longitude / 10000000.0) + ((esekf_double_t)longitudeHp / 1000000000.0);

  int32_t ellipsoid = myGNSS.getElipsoid();
  int8_t ellipsoidHp = myGNSS.getElipsoidHp();
  _CurrentPos.height = ((esekf_double_t)ellipsoid * 10.0 + (esekf_double_t)ellipsoidHp) / 10000.0;
}

void llhToEcef(const PointLLH &llh, PointECEF &ecef)
{
  esekf_double_t latRad = llh.latitude * M_PI_D / 180.0;
  esekf_double_t lonRad = llh.longitude * M_PI_D / 180.0;
  esekf_double_t sinLat = hybrid_math.accurate_sin(latRad);
  esekf_double_t cosLat = hybrid_math.accurate_cos(latRad);

  esekf_double_t N = WGS84_A / sqrt(1.0 - WGS84_E2 * sinLat * sinLat);

  ecef.x = (N + llh.height) * cosLat * hybrid_math.accurate_cos(lonRad);
  ecef.y = (N + llh.height) * cosLat * hybrid_math.accurate_sin(lonRad);
  ecef.z = (N * (1.0 - WGS84_E2) + llh.height) * sinLat;
}

void ecefToEnu_hybrid(const PointECEF &current_ecef, esekf_float_t &east, esekf_float_t &north, esekf_float_t &up)
{
  esekf_double_t dx = current_ecef.x - enu_origin_ecef.x;
  esekf_double_t dy = current_ecef.y - enu_origin_ecef.y;
  esekf_double_t dz = current_ecef.z - enu_origin_ecef.z;

  // Tính toán với double precision, convert sang float ở output
  east = (esekf_float_t)(ecef_to_enu_matrix[0][0] * dx + ecef_to_enu_matrix[0][1] * dy + ecef_to_enu_matrix[0][2] * dz);
  north = (esekf_float_t)(ecef_to_enu_matrix[1][0] * dx + ecef_to_enu_matrix[1][1] * dy + ecef_to_enu_matrix[1][2] * dz);
  up = (esekf_float_t)(ecef_to_enu_matrix[2][0] * dx + ecef_to_enu_matrix[2][1] * dy + ecef_to_enu_matrix[2][2] * dz);
}

void CalculateEnuMatrix_double(const PointLLH &origin_llh)
{
  esekf_double_t latRad = origin_llh.latitude * M_PI_D / 180.0;
  esekf_double_t lonRad = origin_llh.longitude * M_PI_D / 180.0;
  esekf_double_t sLat = hybrid_math.accurate_sin(latRad);
  esekf_double_t cLat = hybrid_math.accurate_cos(latRad);
  esekf_double_t sLon = hybrid_math.accurate_sin(lonRad);
  esekf_double_t cLon = hybrid_math.accurate_cos(lonRad);

  ecef_to_enu_matrix[0][0] = -sLon;
  ecef_to_enu_matrix[0][1] = cLon;
  ecef_to_enu_matrix[0][2] = 0;

  ecef_to_enu_matrix[1][0] = -sLat * cLon;
  ecef_to_enu_matrix[1][1] = -sLat * sLon;
  ecef_to_enu_matrix[1][2] = cLat;

  ecef_to_enu_matrix[2][0] = cLat * cLon;
  ecef_to_enu_matrix[2][1] = cLat * sLon;
  ecef_to_enu_matrix[2][2] = sLat;
}

// --- IMU VARIABLES ---
int i;
esekf_float_t fAcc[3], fGyro[3], fAngle[3];
esekf_float_t accel[2], gyro_z;
esekf_float_t gocIMU = 0.0f;

#define ACC_UPDATE 0x01
#define GYRO_UPDATE 0x02
#define ANGLE_UPDATE 0x04
#define MAG_UPDATE 0x08
#define READ_UPDATE 0x80
static volatile char s_cDataUpdate = 0, s_cCmd = 0xff;

// --- PREDICTION TASK (FLOAT INTERNAL) ---
void predictionTask(void *pvParameters)
{
  const TickType_t xFrequency = 10 / portTICK_PERIOD_MS; // 10 ms => 100Hz
  TickType_t xLastWakeTime = xTaskGetTickCount();
  uint64_t last_time_us = esp_timer_get_time();

  // Pre-allocated arrays (FLOAT)
  static esekf_float_t F[ERROR_STATE_DIM * ERROR_STATE_DIM];
  static esekf_float_t FP[ERROR_STATE_DIM * ERROR_STATE_DIM];
  static esekf_float_t k1[NOMINAL_STATE_DIM], k2[NOMINAL_STATE_DIM],
      k3[NOMINAL_STATE_DIM], k4[NOMINAL_STATE_DIM];

  while (1)
  {
    uint64_t current_time_us = esp_timer_get_time();
    esekf_float_t dt = (current_time_us - last_time_us) / 1000000.0f;
    last_time_us = current_time_us;

    if (dt > 0.05f || dt <= 0)
    {
      vTaskDelayUntil(&xLastWakeTime, xFrequency);
      continue;
    }

    // Đọc IMU
    WitReadReg(AX, 12);
    delay(1);
    while (Serial1.available())
    {
      WitSerialDataIn(Serial1.read());
    }

    esekf_float_t measured_yaw = 0;
    if (s_cDataUpdate)
    {
      for (i = 0; i < 3; i++)
      {
        fAcc[i] = sReg[AX + i] / 32768.0f * 16.0f;
        fGyro[i] = sReg[GX + i] / 32768.0f * 2000.0f;
        fAngle[i] = sReg[Roll + i] / 32768.0f * 180.0f;
      }

      if (s_cDataUpdate & ACC_UPDATE)
      {
        accel[0] = fAcc[0];
        accel[1] = fAcc[1];
        s_cDataUpdate &= ~ACC_UPDATE;
      }
      if (s_cDataUpdate & GYRO_UPDATE)
      {
        gyro_z = fGyro[2] * M_PI_F / 180.0f;
        s_cDataUpdate &= ~GYRO_UPDATE;
      }
      if (s_cDataUpdate & ANGLE_UPDATE)
      {
        measured_yaw = fAngle[2] * M_PI_F / 180.0f;
        gocIMU = fAngle[2];
        s_cDataUpdate &= ~ANGLE_UPDATE;
      }
      s_cDataUpdate = 0;
    }
    else
    {
      vTaskDelayUntil(&xLastWakeTime, xFrequency);
      continue;
    }

    if (xSemaphoreTake(stateMutex, (TickType_t)5) == pdTRUE)
    {
      // =================== PREDICTION (FLOAT) ===================
      esekf_float_t ax_raw = accel[0] - error_state[6];
      esekf_float_t ay_raw = accel[1] - error_state[7];
      esekf_float_t gyro_z_raw = gyro_z - error_state[5];

      esekf_float_t ax_std_body = -ay_raw;
      esekf_float_t ay_std_body = ax_raw;

      esekf_float_t theta = nominal_state[4];
      esekf_float_t cos_th = hybrid_math.fast_cos(theta);
      esekf_float_t sin_th = hybrid_math.fast_sin(theta);

      esekf_float_t forward_accel = ay_std_body;

      // RK4 Integration với float
      k1[0] = nominal_state[2];
      k1[1] = nominal_state[3];
      k1[2] = forward_accel * cos_th;
      k1[3] = forward_accel * sin_th;
      k1[4] = gyro_z_raw;

      esekf_float_t theta_k2 = nominal_state[4] + 0.5f * dt * k1[4];
      k2[0] = nominal_state[2] + 0.5f * dt * k1[2];
      k2[1] = nominal_state[3] + 0.5f * dt * k1[3];
      k2[2] = forward_accel * hybrid_math.fast_cos(theta_k2);
      k2[3] = forward_accel * hybrid_math.fast_sin(theta_k2);
      k2[4] = gyro_z_raw;

      esekf_float_t theta_k3 = nominal_state[4] + 0.5f * dt * k2[4];
      k3[0] = nominal_state[2] + 0.5f * dt * k2[2];
      k3[1] = nominal_state[3] + 0.5f * dt * k2[3];
      k3[2] = forward_accel * hybrid_math.fast_cos(theta_k3);
      k3[3] = forward_accel * hybrid_math.fast_sin(theta_k3);
      k3[4] = gyro_z_raw;

      esekf_float_t theta_k4 = nominal_state[4] + dt * k3[4];
      k4[0] = nominal_state[2] + dt * k3[2];
      k4[1] = nominal_state[3] + dt * k3[3];
      k4[2] = forward_accel * hybrid_math.fast_cos(theta_k4);
      k4[3] = forward_accel * hybrid_math.fast_sin(theta_k4);
      k4[4] = gyro_z_raw;

      for (int i = 0; i < NOMINAL_STATE_DIM; i++)
      {
        nominal_state[i] += (dt / 6.0f) * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]);
      }
      nominal_state[4] = normalizeAngle_float(nominal_state[4]);

      // =================== COVARIANCE PREDICTION (FLOAT) ===================
      memset(F, 0, sizeof(F));
      for (int i = 0; i < ERROR_STATE_DIM; ++i)
        F[i * ERROR_STATE_DIM + i] = 1.0f;

      esekf_float_t theta_pred = nominal_state[4];
      esekf_float_t cos_th_pred = hybrid_math.fast_cos(theta_pred);
      esekf_float_t sin_th_pred = hybrid_math.fast_sin(theta_pred);

      F[0 * ERROR_STATE_DIM + 2] = dt;
      F[1 * ERROR_STATE_DIM + 3] = dt;
      F[2 * ERROR_STATE_DIM + 4] = (-forward_accel * sin_th_pred) * dt;
      F[3 * ERROR_STATE_DIM + 4] = (forward_accel * cos_th_pred) * dt;
      F[4 * ERROR_STATE_DIM + 5] = -dt;
      F[2 * ERROR_STATE_DIM + 6] = -cos_th_pred * dt;
      F[3 * ERROR_STATE_DIM + 6] = -sin_th_pred * dt;

      // Tính FP = F * P (tối ưu hóa với float)
      memset(FP, 0, sizeof(FP));
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        FP[i * ERROR_STATE_DIM + 0] = F[i * ERROR_STATE_DIM + 0] * P[0 * ERROR_STATE_DIM + 0];
        FP[i * ERROR_STATE_DIM + 1] = F[i * ERROR_STATE_DIM + 1] * P[1 * ERROR_STATE_DIM + 1];
        FP[i * ERROR_STATE_DIM + 2] = F[i * ERROR_STATE_DIM + 2] * P[2 * ERROR_STATE_DIM + 2] +
                                      F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 2];
        FP[i * ERROR_STATE_DIM + 3] = F[i * ERROR_STATE_DIM + 3] * P[3 * ERROR_STATE_DIM + 3] +
                                      F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 3];
        FP[i * ERROR_STATE_DIM + 4] = F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 4] +
                                      F[i * ERROR_STATE_DIM + 5] * P[5 * ERROR_STATE_DIM + 4];
        FP[i * ERROR_STATE_DIM + 5] = F[i * ERROR_STATE_DIM + 5] * P[5 * ERROR_STATE_DIM + 5];
        FP[i * ERROR_STATE_DIM + 6] = F[i * ERROR_STATE_DIM + 6] * P[6 * ERROR_STATE_DIM + 6] +
                                      F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 6];
        FP[i * ERROR_STATE_DIM + 7] = F[i * ERROR_STATE_DIM + 7] * P[7 * ERROR_STATE_DIM + 7] +
                                      F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 7];
      }

      // P = F * P * F^T + Q * dt (float)
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < ERROR_STATE_DIM; j++)
        {
          esekf_float_t temp = 0;
          for (int k = 0; k < ERROR_STATE_DIM; k++)
          {
            temp += FP[i * ERROR_STATE_DIM + k] * F[j * ERROR_STATE_DIM + k];
          }
          P[i * ERROR_STATE_DIM + j] = temp + Q[i * ERROR_STATE_DIM + j] * dt;
        }
      }

      // Đảm bảo tính đối xứng
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = i + 1; j < ERROR_STATE_DIM; j++)
        {
          P[j * ERROR_STATE_DIM + i] = P[i * ERROR_STATE_DIM + j];
        }
      }

      xSemaphoreGive(stateMutex);
    }

    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
}

// --- UPDATE TASK (HYBRID: GNSS DOUBLE + ESKF FLOAT) ---
void updateTask(void *pvParameters)
{
  const TickType_t xFrequency = 100 / portTICK_PERIOD_MS;
  TickType_t xLastWakeTime = xTaskGetTickCount();

  // Pre-allocated arrays (FLOAT cho ESKF)
  static esekf_float_t H[5 * ERROR_STATE_DIM] = {0};
  static esekf_float_t R[5 * 5] = {0};
  static esekf_float_t S[5 * 5] = {0};
  static esekf_float_t S_inv[5 * 5];
  static esekf_float_t H_P[5 * ERROR_STATE_DIM] = {0};
  static esekf_float_t P_HT[ERROR_STATE_DIM * 5] = {0};
  static esekf_float_t K[ERROR_STATE_DIM * 5] = {0};
  static esekf_float_t I_minus_KH[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
  static esekf_float_t P_new[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
  static esekf_float_t I_KH[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
  static esekf_float_t P_new_imu[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};

  // Biến tạm cho GNSS (DOUBLE)
  PointECEF tempECEF;

  // Struct để gửi dữ liệu
  ESKF_STATE_DATA eskf_data;

  // Khởi tạo H matrix một lần
  H[0 * ERROR_STATE_DIM + 0] = 1.0f;
  H[1 * ERROR_STATE_DIM + 1] = 1.0f;
  H[2 * ERROR_STATE_DIM + 2] = 1.0f;
  H[3 * ERROR_STATE_DIM + 3] = 1.0f;
  H[4 * ERROR_STATE_DIM + 4] = 1.0f;

  uint8_t rtksta = 0;
  int32_t hAcc = 0, vAcc = 0, headingAcc = 0;

  while (1)
  {
    bool PVT = false;
    esekf_float_t z[5] = {0, 0, 0, 0, 0};

    if (myGNSS.getPVT())
    {
      // Lấy dữ liệu ECEF high precision với DOUBLE
      int32_t ECEFX = myGNSS.getHighResECEFX();
      int8_t ECEFXHp = myGNSS.getHighResECEFXHp();
      int32_t ECEFY = myGNSS.getHighResECEFY();
      int8_t ECEFYHp = myGNSS.getHighResECEFYHp();
      int32_t ECEFZ = myGNSS.getHighResECEFZ();
      int8_t ECEFZHp = myGNSS.getHighResECEFZHp();

      // Tính toán với DOUBLE precision
      tempECEF.x = ((esekf_double_t)ECEFX) / 100.0 + ((esekf_double_t)ECEFXHp) / 10000.0;
      tempECEF.y = ((esekf_double_t)ECEFY) / 100.0 + ((esekf_double_t)ECEFYHp) / 10000.0;
      tempECEF.z = ((esekf_double_t)ECEFZ) / 100.0 + ((esekf_double_t)ECEFZHp) / 10000.0;

      // Chuyển đổi tọa độ với DOUBLE precision
      esekf_float_t current_east, current_north, current_up;
      ecefToEnu_hybrid(tempECEF, current_east, current_north, current_up);

      // Lấy velocity và heading (convert sang float sau khi tính toán)
      esekf_float_t vel_east = ((esekf_float_t)myGNSS.getNedEastVel()) / 1000.0f;
      esekf_float_t vel_north = ((esekf_float_t)myGNSS.getNedNorthVel()) / 1000.0f;
      esekf_float_t heading = ((esekf_float_t)myGNSS.getHeading()) / 100000.0f;

      z[0] = current_east;
      z[1] = current_north;
      z[2] = vel_east;
      z[3] = vel_north;
      z[4] = heading * M_PI_F / 180.0f;
      z[4] = normalizeAngle_float(z[4]);

      rtksta = myGNSS.getCarrierSolutionType();
      hAcc = myGNSS.getHorizontalAccEst();
      vAcc = myGNSS.getVerticalAccEst();
      headingAcc = myGNSS.getHeadingAccEst();

      if (myGNSS.getGroundSpeed() < 300)
        headingAcc = 3500000;

      PVT = true;
    }

    if (PVT && xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
    {
      // =================== MEASUREMENT UPDATE (FLOAT) ===================
      esekf_float_t sigma_p = (esekf_float_t)hAcc / 1000.0f;
      esekf_float_t sigma_v = (esekf_float_t)vAcc / 1000.0f;
      esekf_float_t sigma_yaw = ((esekf_float_t)headingAcc / 10000.0f) * M_PI_F / 180.0f;

      // Adaptive R (FLOAT)
      memset(R, 0, sizeof(R));
      R[0 * 5 + 0] = sigma_p * sigma_p;
      R[1 * 5 + 1] = sigma_p * sigma_p;
      R[2 * 5 + 2] = sigma_v * sigma_v;
      R[3 * 5 + 3] = sigma_v * sigma_v;
      R[4 * 5 + 4] = sigma_yaw * sigma_yaw;

      // Tính S = H * P * H^T + R (FLOAT)
      memset(H_P, 0, sizeof(H_P));
      for (int i = 0; i < 5; i++)
      {
        for (int j = 0; j < ERROR_STATE_DIM; j++)
        {
          for (int k = 0; k < ERROR_STATE_DIM; k++)
          {
            H_P[i * ERROR_STATE_DIM + j] += H[i * ERROR_STATE_DIM + k] * P[k * ERROR_STATE_DIM + j];
          }
        }
      }

      memset(S, 0, sizeof(S));
      for (int i = 0; i < 5; i++)
      {
        for (int j = 0; j < 5; j++)
        {
          for (int k = 0; k < ERROR_STATE_DIM; k++)
          {
            S[i * 5 + j] += H_P[i * ERROR_STATE_DIM + k] * H[j * ERROR_STATE_DIM + k];
          }
          if (i == j)
            S[i * 5 + j] += R[i * 5 + j];
        }
      }

      if (!inverse_5x5_optimized_float(S, S_inv))
      {
        Serial.println("Warning: Matrix inversion failed!");
        xSemaphoreGive(stateMutex);
        vTaskDelayUntil(&xLastWakeTime, xFrequency);
        continue;
      }

      // Tính K = P * H^T * S^-1 (FLOAT)
      memset(P_HT, 0, sizeof(P_HT));
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < 5; j++)
        {
          for (int k = 0; k < ERROR_STATE_DIM; k++)
          {
            P_HT[i * 5 + j] += P[i * ERROR_STATE_DIM + k] * H[j * ERROR_STATE_DIM + k];
          }
        }
      }

      memset(K, 0, sizeof(K));
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < 5; j++)
        {
          for (int k = 0; k < 5; k++)
          {
            K[i * 5 + j] += P_HT[i * 5 + k] * S_inv[k * 5 + j];
          }
        }
      }

      // Innovation và update (FLOAT)
      esekf_float_t innovation[5] = {
          z[0] - nominal_state[0],
          z[1] - nominal_state[1],
          z[2] - nominal_state[2],
          z[3] - nominal_state[3],
          normalizeAngle_float(z[4] - nominal_state[4])};

      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < 5; j++)
        {
          error_state[i] += K[i * 5 + j] * innovation[j];
        }
      }

      // Inject error và reset (FLOAT)
      nominal_state[0] += error_state[0];
      nominal_state[1] += error_state[1];
      nominal_state[2] += error_state[2];
      nominal_state[3] += error_state[3];
      nominal_state[4] += error_state[4];
      nominal_state[4] = normalizeAngle_float(nominal_state[4]);

      memset(error_state, 0, 4 * sizeof(esekf_float_t));

      // Update covariance: P = (I - K*H) * P (FLOAT)
      memset(I_minus_KH, 0, sizeof(I_minus_KH));
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < ERROR_STATE_DIM; j++)
        {
          I_minus_KH[i * ERROR_STATE_DIM + j] = (i == j) ? 1.0f : 0.0f;
          for (int k = 0; k < 5; k++)
          {
            I_minus_KH[i * ERROR_STATE_DIM + j] -= K[i * 5 + k] * H[k * ERROR_STATE_DIM + j];
          }
        }
      }

      memset(P_new, 0, sizeof(P_new));
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < ERROR_STATE_DIM; j++)
        {
          for (int k = 0; k < ERROR_STATE_DIM; k++)
          {
            P_new[i * ERROR_STATE_DIM + j] += I_minus_KH[i * ERROR_STATE_DIM + k] * P[k * ERROR_STATE_DIM + j];
          }
        }
      }
      memcpy(P, P_new, sizeof(P_new));

      // Đảm bảo tính đối xứng
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = i + 1; j < ERROR_STATE_DIM; j++)
        {
          P[j * ERROR_STATE_DIM + i] = P[i * ERROR_STATE_DIM + j];
        }
      }

      // =================== BƯỚC 4: MINI-UPDATE DÙNG IMU YAW (RELATIVE ANGLE) ===================
      // Mục đích: Dùng độ mượt của IMU Yaw để hiệu chỉnh độ trôi (b_gz) và làm mượt dpsi
      if (fabsf(C0_yaw) > 1e-4f) // Chỉ thực hiện nếu C0_yaw đã được tính toán ở setup
      {
        // 1. Góc IMU đã được chuyển sang radian (gocIMU đang là độ)
        esekf_float_t measured_imu_yaw = gocIMU * M_PI_F / 180.0f;

        // 2. Áp dụng góc bù trừ C0_yaw
        // Đây là góc IMU đã được chuyển về hệ Navigation (North-referenced)
        esekf_float_t corrected_imu_yaw = normalizeAngle_float(measured_imu_yaw - C0_yaw);

        // 3. Tính toán Innovation (Sai số đo lường)
        // Sai số = corrected_imu_yaw - nominal_state[4] (yaw hiện tại)
        esekf_float_t innovation_imu_yaw = normalizeAngle_float(corrected_imu_yaw - nominal_state[4]);

        // Ma trận H cho phép đo yaw (1x8)
        // H_yaw = [0, 0, 0, 0, 1, 0, 0, 0]
        // Ma trận R cho phép đo yaw (1x1) là R_yaw (đã được tính ở setup, R_yaw nhỏ)

        // S_yaw = H_yaw * P * H_yaw^T + R_yaw  (S là số vô hướng 1x1)
        esekf_float_t S_yaw = P[4 * ERROR_STATE_DIM + 4] + R_yaw;

        // K_yaw = P * H_yaw^T * S_yaw_inv (K là véc-tơ 8x1)
        esekf_float_t K_yaw[ERROR_STATE_DIM];
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          // P[i * ERROR_STATE_DIM + 4] là cột 4 của P, tương đương P * H^T
          K_yaw[i] = P[i * ERROR_STATE_DIM + 4] / S_yaw;
        }

        // Cập nhật trạng thái lỗi (cộng dồn với kết quả từ GPS update)
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          error_state[i] += K_yaw[i] * innovation_imu_yaw;
        }

        // Cập nhật hiệp phương sai: P = (I - K * H) * P (Sequential Update)
        memset(I_KH, 0, sizeof(I_KH));
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          for (int j = 0; j < ERROR_STATE_DIM; j++)
          {
            // H_yaw chỉ có phần tử thứ 4 là khác 0 (trạng thái dpsi)
            esekf_float_t kh_ij = K_yaw[i] * ((j == 4) ? 1.0f : 0.0f);
            I_KH[i * ERROR_STATE_DIM + j] = ((i == j) ? 1.0f : 0.0f) - kh_ij;
          }
        }

        memset(P_new_imu, 0, sizeof(P_new_imu));
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          for (int j = 0; j < ERROR_STATE_DIM; j++)
          {
            for (int k = 0; k < ERROR_STATE_DIM; k++)
            {
              P_new_imu[i * ERROR_STATE_DIM + j] += I_KH[i * ERROR_STATE_DIM + k] * P[k * ERROR_STATE_DIM + j];
            }
          }
        }
        memcpy(P, P_new_imu, sizeof(P)); // Ghi đè P mới

        // Đảm bảo tính đối xứng
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          for (int j = i + 1; j < ERROR_STATE_DIM; j++)
          {
            P[j * ERROR_STATE_DIM + i] = P[i * ERROR_STATE_DIM + j];
          }
        }

        // Debug IMU yaw update
        Serial.print("IMU_YAW_UPDATE: corrected=");
        Serial.print(corrected_imu_yaw * 180.0f / M_PI_F, 2);
        Serial.print("°, innovation=");
        Serial.print(innovation_imu_yaw * 180.0f / M_PI_F, 2);
        Serial.print("°, C0=");
        Serial.print(C0_yaw * 180.0f / M_PI_F, 2);
        Serial.println("°");
      } // Kết thúc Update IMU

      // =================== GỬI DỮ LIỆU QUA SERIAL ===================
      // Điền dữ liệu vào struct
      eskf_data.timestamp = millis();

      // Copy nominal state
      memcpy(eskf_data.state, nominal_state, sizeof(nominal_state));

      // Copy đường chéo của ma trận P
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        eskf_data.P_diag[i] = P[i * ERROR_STATE_DIM + i];
      }

      // Dữ liệu cảm biến thô
      eskf_data.raw_data[0] = accel[0];                           // ax
      eskf_data.raw_data[1] = accel[1];                           // ay
      eskf_data.raw_data[2] = gocIMU;                             // yaw IMU (độ)
      eskf_data.raw_data[3] = gyro_z;                             // gz (rad/s)
      eskf_data.raw_data[4] = (esekf_float_t)tempECEF.x;          // ECEF X
      eskf_data.raw_data[5] = (esekf_float_t)tempECEF.y;          // ECEF Y
      eskf_data.raw_data[6] = nominal_state[2];                   // v_east
      eskf_data.raw_data[7] = nominal_state[3];                   // v_north
      eskf_data.raw_data[8] = nominal_state[4] * 180.0f / M_PI_F; // heading (độ)

      // Độ chính xác GNSS
      eskf_data.gnss_accuracy[0] = (esekf_float_t)hAcc / 1000.0f;        // hAcc (m)
      eskf_data.gnss_accuracy[1] = (esekf_float_t)vAcc / 1000.0f;        // vAcc (m/s)
      eskf_data.gnss_accuracy[2] = (esekf_float_t)headingAcc / 10000.0f; // headingAcc (độ)

      // Gửi qua Serial dưới dạng binary
      Serial.write(&SYNC_BYTE, 1); // Sync 0xAA
      Serial.write((uint8_t *)&eskf_data, sizeof(ESKF_STATE_DATA));
      Serial.flush();

      // DEBUG: In ra Serial dưới dạng text để kiểm tra
      // Serial.print("ESKF_DATA,");
      // Serial.print(eskf_data.timestamp);
      // Serial.print(",");
      // for (int i = 0; i < 5; i++)
      // {
      //   Serial.print(eskf_data.state[i], 6);
      //   Serial.print(",");
      // }
      // for (int i = 0; i < 8; i++)
      // {
      //   Serial.print(eskf_data.P_diag[i], 8);
      //   Serial.print(",");
      // }
      // for (int i = 0; i < 9; i++)
      // {
      //   Serial.print(eskf_data.raw_data[i], 4);
      //   Serial.print(",");
      // }
      // for (int i = 0; i < 3; i++)
      // {
      //   Serial.print(eskf_data.gnss_accuracy[i], 4);
      //   if (i < 2)
      //     Serial.print(",");
      // }
      // Serial.println();

      xSemaphoreGive(stateMutex);
    }

    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
}

// --- IMU COMMUNICATION FUNCTIONS ---
static void SensorUartSend(uint8_t *p_data, uint32_t uiSize)
{
  Serial1.write(p_data, uiSize);
  Serial1.flush();
}

static void CopeSensorData(uint32_t uiReg, uint32_t uiRegNum)
{
  for (int i = 0; i < uiRegNum; i++)
  {
    switch (uiReg)
    {
    case AZ:
      s_cDataUpdate |= ACC_UPDATE;
      break;
    case GZ:
      s_cDataUpdate |= GYRO_UPDATE;
      break;
    case HZ:
      s_cDataUpdate |= MAG_UPDATE;
      break;
    case Yaw:
      s_cDataUpdate |= ANGLE_UPDATE;
      break;
    default:
      s_cDataUpdate |= READ_UPDATE;
      break;
    }
    uiReg++;
  }
}

static void Delayms(uint16_t ucMs)
{
  delay(ucMs);
}

void setup()
{
  // KHỐI 1: KHỞI TẠO HỆ THỐNG CỐT LÕI
  setCpuFrequencyMhz(240);
  delay(500); // Chờ hệ thống ổn định

  Serial.begin(115200);
  Serial.println("System Core Initializing...");

  // Tạo Mutex để bảo vệ các biến dùng chung giữa các task
  stateMutex = xSemaphoreCreateMutex();
  if (stateMutex == NULL)
  {
    Serial.println("Error: Failed to create mutex!");
    while (1)
      ; // Treo hệ thống nếu không tạo được Mutex
  }

  // KHỐI 2: KHỞI TẠO CÁC GIAO TIẾP PHẦN CỨNG
  Serial.println("Initializing Hardware Interfaces...");

  // Giao tiếp UART với IMU
  Serial1.begin(115200, SERIAL_8N1, RXD2, TXD2);

  // Giao tiếp I2C với module GNSS
  Wire.begin(21, 22);
  delay(500); // Chờ GNSS khởi động

  // KHỐI 3: KHỞI TẠO CÁC CẢM BIẾN
  Serial.println("Initializing Sensors...");

  // Khởi tạo IMU WitMotion
  WitInit(WIT_PROTOCOL_MODBUS, 0x50);
  WitSerialWriteRegister(SensorUartSend);
  WitRegisterCallBack(CopeSensorData);
  WitDelayMsRegister(Delayms);
  // WitSetHeadingtoZero();

  // Khởi tạo Module GNSS u-blox
  if (myGNSS.begin(Wire, 0x42) == false)
  {
    Serial.println(F("u-blox GNSS not detected. Freezing."));
    while (1)
      ;
  }
  myGNSS.setI2COutput(COM_TYPE_UBX);
  myGNSS.setNavigationFrequency(10); // 5Hz

  // KHỐI 4: CHỜ TÍN HIỆU VÀ THIẾT LẬP GỐC TỌA ĐỘ ENU
  Serial.println("Waiting for RTK Fixed solution to set ENU origin...");
  bool origin_set = false;
  while (!origin_set)
  {
    uint8_t fixType = 0;
    uint8_t carrierSolution = 0;

    if (myGNSS.getPVT())
    {
      fixType = myGNSS.getFixType();
      carrierSolution = myGNSS.getCarrierSolutionType();

      Serial.print("Fix Type: ");
      Serial.print(fixType);
      Serial.print(", Carrier Solution: ");
      Serial.println(carrierSolution);

      // Khi có tín hiệu RTK Fixed, thiết lập gốc tọa độ
      if (fixType >= 3 && carrierSolution >= 1)
      {
        Serial.println("RTK Fixed acquired. Setting ENU origin...");

        // Lấy tọa độ gốc LLH từ GNSS để tính ma trận xoay
        // GetRTK(initialLLH);
        initialLLH.latitude = 22.626686000;
        initialLLH.longitude = 120.266022880;
        initialLLH.height = 24.2480;
        CalculateEnuMatrix_double(initialLLH);

        // Lấy tọa độ gốc ECEF trực tiếp từ GNSS
        enu_origin_ecef.x = ((esekf_double_t)myGNSS.getHighResECEFX() / 100.0) + ((esekf_double_t)myGNSS.getHighResECEFXHp() / 10000.0);
        enu_origin_ecef.y = ((esekf_double_t)myGNSS.getHighResECEFY() / 100.0) + ((esekf_double_t)myGNSS.getHighResECEFYHp() / 10000.0);
        enu_origin_ecef.z = ((esekf_double_t)myGNSS.getHighResECEFZ() / 100.0) + ((esekf_double_t)myGNSS.getHighResECEFZHp() / 10000.0);

        Serial.print("ENU Origin ECEF: ");
        Serial.print(enu_origin_ecef.x, 4);
        Serial.print(", ");
        Serial.print(enu_origin_ecef.y, 4);
        Serial.print(", ");
        Serial.println(enu_origin_ecef.z, 4);

        origin_set = true;
      }
    }
    delay(200); // Chờ 200ms trước khi thử lại
  }

  // KHỐI 5: KHỞI TẠO TRẠNG THÁI VÀ HIỆP PHƯƠNG SAI CHO ES-EKF
  Serial.println("Initializing ES-EKF State and Covariance...");
  if (myGNSS.getNAVHPPOSECEF())
  {
    uint32_t hAcc_p = myGNSS.getHorizontalAccEst(); // mm
    uint32_t vAcc_p = myGNSS.getVerticalAccEst();   // mm/s (giả định)
    uint32_t headingAcc_p;                          // độ (giả định)                                                         // 1e-5 độ

    esekf_float_t sigma_p = (esekf_float_t)hAcc_p / 1000.0f;                     // mm -> m
    esekf_float_t sigma_v = vAcc_p > 0 ? (esekf_float_t)vAcc_p / 1000.0f : 0.1f; // mm/s -> m/s, mặc định 0.1 m/s nếu không có
    esekf_float_t sigma_yaw;                                                     // độ -> radian, mặc định 5°
    esekf_float_t sigma_bg_z = 0.1f;                                             // rad/s (bias con quay)
    esekf_float_t sigma_ba = 0.3f;
    bool heading_valid = false;
    while (!heading_valid)
    {
      if (myGNSS.getGroundSpeed() > V_MIN_HEADING * 1000)
      {
        heading_valid = true;
        headingAcc_p = myGNSS.getHeadingAccEst();
        sigma_yaw = headingAcc_p > 0 ? (esekf_float_t)headingAcc_p / 10000.0f * M_PI_F / 180.0f : 0.0873f; // độ -> radian, mặc định 5°
        // sigma_yaw = 5.0 * M_PI / 180.0; // Đặt độ chính xác hướng rất lớn (5 độ)
        nominal_state[4] = myGNSS.getHeading() / 100000.0f * M_PI_F / 180.0f;
        nominal_state[4] = normalizeAngle_float(nominal_state[4]);
      }
      else
      {
        Serial.println("Waiting for valid GNSS heading (ground speed too low)...");
        delay(50);
      }
    }
    Serial.print("Initial Heading from GNSS (rad): ");
    Serial.println(nominal_state[4]);

    Serial.print("Initial Position Std Dev (m): ");
    Serial.print(sigma_p, 4);
    Serial.print(", Velocity Std Dev (m/s): ");
    Serial.print(sigma_v, 4);
    Serial.print(", Heading Std Dev (deg): ");
    Serial.print(sigma_yaw * 180.0f / M_PI_F, 4);

    P[0 * ERROR_STATE_DIM + 0] = sigma_p * sigma_p;       // δp_E
    P[1 * ERROR_STATE_DIM + 1] = sigma_p * sigma_p;       // δp_N
    P[2 * ERROR_STATE_DIM + 2] = sigma_v * sigma_v;       // δv_E
    P[3 * ERROR_STATE_DIM + 3] = sigma_v * sigma_v;       // δv_N
    P[4 * ERROR_STATE_DIM + 4] = sigma_yaw * sigma_yaw;   // δψ
    P[5 * ERROR_STATE_DIM + 5] = sigma_bg_z * sigma_bg_z; // b_g,z
    P[6 * ERROR_STATE_DIM + 6] = sigma_ba * sigma_ba;     // b_a,x
    P[7 * ERROR_STATE_DIM + 7] = sigma_ba * sigma_ba;     // b_a,y
  }

  // Khởi tạo trạng thái danh nghĩa (Nominal State) ---
  nominal_state[0] = 0.0; // Vị trí East ban đầu tại gốc
  nominal_state[1] = 0.0; // Vị trí North ban đầu tại gốc
  nominal_state[2] = 0.0; // Vận tốc East ban đầu
  nominal_state[3] = 0.0; // Vận tốc North ban đầu

  // Lấy heading từ IMU và tính toán góc lệch C0_yaw ---
  WitReadReg(AX, 12);
  delay(10);
  while (Serial1.available())
  {
    WitSerialDataIn(Serial1.read());
  }

  if (s_cDataUpdate & ANGLE_UPDATE)
  {
    esekf_float_t initial_imu_yaw_deg = sReg[Yaw] / 32768.0f * 180.0f;
    C0_yaw = initial_imu_yaw_deg * M_PI_F / 180.0f - nominal_state[4];
    s_cDataUpdate &= ~ANGLE_UPDATE;
    Serial.print("Initial IMU yaw offset (C0_yaw) calculated: ");
    Serial.println(C0_yaw * 180.0f / M_PI_F);
  }
  else
  {
    C0_yaw = 0.0f;
    Serial.println("Warning: Could not read initial IMU heading.");
  }

  // --- 5.4: Khởi tạo ma trận nhiễu quá trình Q ---
  esekf_float_t vel_noise_std_dev = 1.5f;
  esekf_float_t accel_bias_noise_std_dev = 0.86999f;
  esekf_float_t gyro_bias_noise_std_dev = 0.8994f;

  Q[2 * ERROR_STATE_DIM + 2] = vel_noise_std_dev * vel_noise_std_dev;
  Q[3 * ERROR_STATE_DIM + 3] = vel_noise_std_dev * vel_noise_std_dev;
  Q[4 * ERROR_STATE_DIM + 4] = (100.0f * M_PI_F / 180.0f) * (100.0f * M_PI_F / 180.0f); // Nhiễu đo của Gyro
  Q[5 * ERROR_STATE_DIM + 5] = gyro_bias_noise_std_dev * gyro_bias_noise_std_dev;
  Q[6 * ERROR_STATE_DIM + 6] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;
  Q[7 * ERROR_STATE_DIM + 7] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;

  // Khởi tạo ma trận nhiễu đo lường R cho IMU Yaw ---
  esekf_float_t yaw_measurement_std_dev_rad = 0.5f * M_PI_F / 180.0f; // Sai số đo của IMU yaw là 1.1 độ
  R_yaw = yaw_measurement_std_dev_rad * yaw_measurement_std_dev_rad;

  // KHỐI 6: TẠO CÁC TÁC VỤ (TASKS)
  Serial.println("Creating FreeRTOS tasks...");
  xTaskCreatePinnedToCore(predictionTask, "Prediction", 8192, NULL, 5, NULL, 0); // Core 0
  xTaskCreatePinnedToCore(updateTask, "Update", 12288, NULL, 5, NULL, 1);        // Core 1
  Serial.println("System initialized successfully. ES-EKF is running.");
  delay(2000);
}

void loop()
{
  vTaskDelay(1000 / portTICK_PERIOD_MS);
}