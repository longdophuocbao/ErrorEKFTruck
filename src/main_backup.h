// ES-EKF
#include <Wire.h>
#include <SparkFun_u-blox_GNSS_v3.h>
// #include <ModbusMaster.h>
// #include <WTGAHRS3_485.h>
#include <esp_now.h>
#include <WiFi.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <REG.h>
#include <wit_c_sdk.h>

// #include "BluetoothSerial.h"
#include <freertos/semphr.h> // Để dùng Mutex

// BluetoothSerial SerialBT;

// --- CẤU HÌNH BỘ LỌC ---
const int NOMINAL_STATE_DIM = 5; // Kích thước của Nominal State
const int ERROR_STATE_DIM = 8;   // Kích thước của Error State

const float V_MIN_HEADING = 1.0;

struct ESKF_STATE_DATA
{
  // Vector trạng thái Nominal: 5 phần tử
  float state_x[NOMINAL_STATE_DIM];

  // Đường chéo Ma trận Hiệp phương sai P: 8 phần tử
  float P_diag[ERROR_STATE_DIM];

  // Dữ liệu cảm biến thô (để debug, 9 phần tử: ax, ay, yawIMU, gz, ecef_x, ecef_y, vn, ve, heading)
  float raw_data[9];

  // Dữ liệu độ chính xác GNSS (để debug, 3 phần tử: hAcc, vAcc, sAcc)
  float gnss_accuracy[3];

  uint32_t timestamp;
} __attribute__((packed)); // Truyền binary tối ưu

// Cấu trúc để nhận THAM SỐ TINH CHỈNH (Giữ nguyên, ví dụ 3xQ, 3xR)
struct TUNING_CONFIG
{
  float Q_noise_params[3];
  float R_noise_params[3];
} __attribute__((packed));

TUNING_CONFIG currentConfig;

#define RXD2 16
#define TXD2 17
#define SENSOR_SLAVE_ID 0x50

int i;
float fAcc[3], fGyro[3], fAngle[3];

#define ACC_UPDATE 0x01
#define GYRO_UPDATE 0x02
#define ANGLE_UPDATE 0x04
#define MAG_UPDATE 0x08
#define READ_UPDATE 0x80
static volatile char s_cDataUpdate = 0, s_cCmd = 0xff;

// --- BIẾN TOÀN CỤC ---
SFE_UBLOX_GNSS myGNSS;

// ModbusMaster node;
// HardwareSerial RS485(1);
// WTGAHRS3_485 sensor(node);
SemaphoreHandle_t stateMutex;

// --- KHAI BÁO CÁC STRUCT VÀ HÀM TIỆN ÍCH (Giữ nguyên) ---
struct PointLLH
{
  double latitude, longitude, height;
};
struct PointECEF
{
  double x, y, z;
};

PointLLH initialLLH;
PointECEF initialECEF;

// Trạng thái danh nghĩa (Nominal State) - Ước tính chính của hệ thống
double nominal_state[NOMINAL_STATE_DIM] = {0}; // [p_e, p_n, v_e, v_n, theta]
double error_state[ERROR_STATE_DIM] = {0};     // [dp_e, dp_n, dv_e, dv_n, dtheta, db_gz, db_ax, db_ay]
// Ma trận hiệp phương sai của TRẠNG THÁI LỖI (Error State Covariance)
double P[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};

// Ma trận nhiễu quá trình cho TRẠNG THÁI LỖI
double Q[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};

double R_yaw;
double C0_yaw = 0.0;
// --- Các biến phần cứng và hệ quy chiếu ---
double accel[2], gyro_z;
PointECEF enu_origin_ecef;
double ecef_to_enu_matrix[3][3];

const double WGS84_A = 6378137.0;
const double WGS84_E2 = 0.00669437999014;

/**
 * @brief Lấy dữ liệu RTK từ module GNSS và lưu vào cấu trúc PointLLH.
 * @param _CurrentPos Tham chiếu đến cấu trúc PointLLH để lưu dữ liệu RTK.
 */
void GetRTK(PointLLH &_CurrentPos)
{
  // put your main code here, to run repeatedly:
  int32_t latitude = myGNSS.getHighResLatitude();
  int8_t latitudeHp = myGNSS.getHighResLatitudeHp();
  _CurrentPos.latitude = ((double)latitude / 10000000.0) + ((double)latitudeHp / 1000000000.0);

  int32_t longitude = myGNSS.getHighResLongitude();
  int8_t longitudeHp = myGNSS.getHighResLongitudeHp();
  _CurrentPos.longitude = ((double)longitude / 10000000.0) + ((double)longitudeHp / 1000000000.0);

  int32_t ellipsoid = myGNSS.getElipsoid();
  int8_t ellipsoidHp = myGNSS.getElipsoidHp();
  _CurrentPos.height = ((double)ellipsoid * 10.0 + (double)ellipsoidHp) / 10000.0;
}

/**
 * @brief Chuyển đổi tọa độ từ LLH (Vĩ độ, Kinh độ, Cao độ) sang ECEF (X, Y, Z).
 * @param llh Điểm đầu vào ở dạng LLH (đơn vị: độ, độ, mét).
 * @param ecef Điểm đầu ra sẽ được lưu vào đây.
 */
void llhToEcef(const PointLLH &llh, PointECEF &ecef)
{
  double latRad = llh.latitude * M_PI / 180.0;
  double lonRad = llh.longitude * M_PI / 180.0;
  double sinLat = sin(latRad);
  double cosLat = cos(latRad);

  double N = WGS84_A / sqrt(1.0 - WGS84_E2 * sinLat * sinLat);

  ecef.x = (N + llh.height) * cosLat * cos(lonRad);
  ecef.y = (N + llh.height) * cosLat * sin(lonRad);
  ecef.z = (N * (1.0 - WGS84_E2) + llh.height) * sinLat;
}

/**
 * @brief Thiết lập ma trận chuyển đổi từ ECEF sang ENU dựa trên gốc ENU.
 * @param origin_llh Tọa độ gốc ENU ở dạng LLH.
 * @param origin_ecef Tọa độ gốc ENU ở dạng ECEF.
 */
void ecefToEnu(const PointECEF &current_ecef, double &east, double &north, double &up)
{
  double dx = current_ecef.x - enu_origin_ecef.x;
  double dy = current_ecef.y - enu_origin_ecef.y;
  double dz = current_ecef.z - enu_origin_ecef.z;
  // Serial.print("dx=");
  // Serial.print(dx, 4);
  // Serial.print(",dy=");
  // Serial.println(dy, 4);

  east = ecef_to_enu_matrix[0][0] * dx + ecef_to_enu_matrix[0][1] * dy + ecef_to_enu_matrix[0][2] * dz;
  north = ecef_to_enu_matrix[1][0] * dx + ecef_to_enu_matrix[1][1] * dy + ecef_to_enu_matrix[1][2] * dz;
  up = ecef_to_enu_matrix[2][0] * dx + ecef_to_enu_matrix[2][1] * dy + ecef_to_enu_matrix[2][2] * dz;
}

/**
 * @brief Tính ma trận chuyển đổi từ ECEF sang ENU dựa trên gốc ENU.
 * @param origin_llh Tọa độ gốc ENU ở dạng LLH.
 */
void CalculateEnuMatrix(const PointLLH &origin_llh)
{
  double latRad = origin_llh.latitude * M_PI / 180.0;
  double lonRad = origin_llh.longitude * M_PI / 180.0;
  double sLat = sin(latRad);
  double cLat = cos(latRad);
  double sLon = sin(lonRad);
  double cLon = cos(lonRad);

  // Ma trận xoay từ ECEF sang ENU
  ecef_to_enu_matrix[0][0] = -sLon;
  ecef_to_enu_matrix[0][1] = cLon;
  ecef_to_enu_matrix[0][2] = 0;

  ecef_to_enu_matrix[1][0] = -sLat * cLon;
  ecef_to_enu_matrix[1][1] = -sLat * sLon;
  ecef_to_enu_matrix[1][2] = cLat;

  ecef_to_enu_matrix[2][0] = cLat * cLon;
  ecef_to_enu_matrix[2][1] = cLat * sLon;
  ecef_to_enu_matrix[2][2] = sLat;

  Serial.println("ECEF to ENU matrix:");
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      Serial.print(ecef_to_enu_matrix[i][j], 6);
      Serial.print(" ");
    }
    Serial.println();
  }
}

/**
 * @brief Chuẩn hóa góc về khoảng [-π, π].
 * @param angle Góc đầu vào (radian).
 * @return Góc đã chuẩn hóa (radian).
 */
double normalizeAngle(double angle)
{
  return atan2(sin(angle), cos(angle));
}

// --- KẾT THÚC KHAI BÁO ---
const int LOOP_FREQ_HZ = 90;                              // Tần số vòng lặp mong muốn (Hz)
const unsigned long LOOP_PERIOD_MS = 1000 / LOOP_FREQ_HZ; // Chu kỳ vòng lặp (ms)
void loopRate(int freq)
{
  static unsigned long lastLoopTime = 0;
  unsigned long currentTime = millis();
  unsigned long loopPeriodMs = 1000 / freq; // Chu kỳ vòng lặp (ms)

  // Tính thời gian đã trôi qua kể từ lần chạy trước
  unsigned long elapsedTime = currentTime - lastLoopTime;

  // Nếu thời gian đã trôi qua nhỏ hơn chu kỳ mong muốn, thêm độ trễ
  if (elapsedTime < loopPeriodMs)
  {
    // Serial.println("adfdasfd");
    TickType_t delayTicks = (loopPeriodMs - elapsedTime) / portTICK_PERIOD_MS;
    vTaskDelay(delayTicks);
  }

  // Cập nhật thời gian cho lần chạy tiếp theo
  lastLoopTime = millis();
}

// Sử dụng LU decomposition thay vì Cholesky cho ổn định hơn
bool inverse_5x5_optimized(const double *A, double *A_inv)
{
  double LU[25];
  memcpy(LU, A, sizeof(LU));

  int pivot[5];

  // LU decomposition
  for (int i = 0; i < 5; i++)
  {
    // Tìm pivot
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
      // Swap rows
      for (int j = 0; j < 5; j++)
      {
        double temp = LU[i * 5 + j];
        LU[i * 5 + j] = LU[max_row * 5 + j];
        LU[max_row * 5 + j] = temp;
      }
    }

    pivot[i] = max_row;

    // Elimination
    for (int j = i + 1; j < 5; j++)
    {
      LU[j * 5 + i] /= LU[i * 5 + i];
      for (int k = i + 1; k < 5; k++)
      {
        LU[j * 5 + k] -= LU[j * 5 + i] * LU[i * 5 + k];
      }
    }
  }

  // Forward/back substitution cho từng cột
  for (int col = 0; col < 5; col++)
  {
    double x[5] = {0};

    // Forward
    for (int i = 0; i < 5; i++)
    {
      x[i] = (col == pivot[i]) ? 1.0f : 0.0f;
      for (int j = 0; j < i; j++)
      {
        x[i] -= LU[i * 5 + j] * x[j];
      }
    }

    // Back
    for (int i = 4; i >= 0; i--)
    {
      for (int j = i + 1; j < 5; j++)
      {
        x[i] -= LU[i * 5 + j] * x[j];
      }
      x[i] /= LU[i * 5 + i];
    }

    // Copy kết quả
    for (int i = 0; i < 5; i++)
    {
      A_inv[i * 5 + col] = x[i];
    }
  }

  return true;
}

bool inverse_5x5(const double *A, double *A_inv)
{
  const int N = 5;
  double L[N][N] = {0};
  double y[N], x[N];
  double I[N][N] = {0};

  // Copy A vào mảng 2D để dễ xử lý
  double M[N][N];
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      M[i][j] = A[i * N + j];

  // --- Bước 1: Cholesky decomposition: M = L * L^T ---
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      double sum = M[i][j];
      for (int k = 0; k < j; k++)
        sum -= L[i][k] * L[j][k];

      if (i == j)
      {
        if (sum <= 1e-12)
          return false; // không SPD
        L[i][j] = sqrt(sum);
      }
      else
      {
        L[i][j] = sum / L[j][j];
      }
    }
  }

  // --- Bước 2: Giải hệ cho từng cột của ma trận nghịch đảo ---
  for (int col = 0; col < N; col++)
  {
    // vector đơn vị e_col
    for (int i = 0; i < N; i++)
      I[i][col] = (i == col) ? 1.0 : 0.0;

    // Ly = e_col (forward substitution)
    for (int i = 0; i < N; i++)
    {
      double sum = I[i][col];
      for (int k = 0; k < i; k++)
        sum -= L[i][k] * y[k];
      y[i] = sum / L[i][i];
    }

    // L^T x = y (back substitution)
    for (int i = N - 1; i >= 0; i--)
    {
      double sum = y[i];
      for (int k = i + 1; k < N; k++)
        sum -= L[k][i] * x[k];
      x[i] = sum / L[i][i];
    }

    // Gán cột kết quả vào A_inv
    for (int i = 0; i < N; i++)
      A_inv[i * N + col] = x[i];
  }

  return true;
}
float gocIMU = 0.0;
void predictionTask(void *pvParameters)
{
  // Tần số mục tiêu là 98Hz -> Chu kỳ = 1000/98 ≈ 10.2ms.
  // Làm tròn đến số nguyên tick gần nhất (1 tick = 1ms) là 10ms -> Tần số thực tế 100Hz.
  const TickType_t xFrequency = 10 / portTICK_PERIOD_MS;
  TickType_t xLastWakeTime = xTaskGetTickCount();

  // (Đề xuất) Đo lường dt thực tế để tăng độ chính xác
  uint64_t last_time_us = esp_timer_get_time();
  long prediction_counter = 0;
  while (1)
  {
    // Tính toán dt thực tế giữa các lần lặp
    uint64_t current_time_us = esp_timer_get_time();
    double dt = (current_time_us - last_time_us) / 1000000.0;
    last_time_us = current_time_us;
    double measured_yaw = 0;
    // Bỏ qua nếu dt quá lớn (ví dụ: do một tác vụ khác gây trễ)
    if (dt > 0.05 || dt <= 0)
    {
      vTaskDelayUntil(&xLastWakeTime, xFrequency); // Vẫn phải delay để tránh task chạy liên tục
      continue;
    }
    WitReadReg(AX, 12);
    delay(10);
    while (Serial1.available())
    {
      WitSerialDataIn(Serial1.read());
    }
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
        gyro_z = fGyro[2] * M_PI / 180.0; // Chuyển sang rad/s
        s_cDataUpdate &= ~GYRO_UPDATE;
      }
      if (s_cDataUpdate & ANGLE_UPDATE)
      {
        // Serial.print("angle:");
        // Serial.print(fAngle[0], 3);
        // Serial.print(" ");
        // Serial.print(fAngle[1], 3);
        // Serial.print(" ");
        measured_yaw = fAngle[2] * M_PI / 180.0;
        // Serial.print("\r\n");
        s_cDataUpdate &= ~ANGLE_UPDATE;
      }

      s_cDataUpdate = 0;
    }
    else
    {
      // Nếu không có dữ liệu mới, bỏ qua lần lặp này
      vTaskDelayUntil(&xLastWakeTime, xFrequency);
      continue;
    }
    /// Hiển thị dữ liệu cảm biến đọc được (debug)
    /*
    Serial.print("dt=");
    Serial.print(dt, 4);
    Serial.print(", ax=");
    Serial.print(accel[0], 3);
    Serial.print(", ay=");
    Serial.print(accel[1], 3);
    Serial.print(", gz=");
    Serial.print(gyro_z, 3);
    Serial.print(", yaw=");
    Serial.println(measured_yaw*180.0/M_PI, 3);*/

    if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
    {
      // =================== BƯỚC 1: CẬP NHẬT TRẠNG THÁI DANH NGHĨA (NOMINAL STATE) ===================
      // Sử dụng RK4 để tích hợp IMU vào nominal_state, làm cho nó bị trôi đi một cách tự nhiên.

      double ax_raw = accel[0] - error_state[6];
      double ay_raw = accel[1] - error_state[7];
      double gyro_z_raw = gyro_z - error_state[5];
      gocIMU = fAngle[2];

      // Chuyển gia tốc từ hệ toạ độ body sang hệ body tiêu chuẩn ( y trước, x phải)
      double ax_std_body = -ay_raw;
      double ay_std_body = ax_raw;

      double theta = nominal_state[4];
      double cos_th = cos(theta);
      double sin_th = sin(theta);

      double forward_accel = ay_std_body;

      // Serial.print(forward_accel, 4);
      // forward_accel = 0.0;
      double ax_nav_east = forward_accel * cos_th;
      double ay_nav_north = forward_accel * sin_th;

      double k1[NOMINAL_STATE_DIM], k2[NOMINAL_STATE_DIM], k3[NOMINAL_STATE_DIM], k4[NOMINAL_STATE_DIM];

      // -- Tính k1
      k1[0] = nominal_state[2];
      k1[1] = nominal_state[3];
      k1[2] = ax_nav_east;
      k1[3] = ay_nav_north;
      k1[4] = gyro_z_raw;

      // --Tính k2
      double theta_k2 = nominal_state[4] + 0.5 * dt * k1[4];
      k2[0] = nominal_state[2] + 0.5 * dt * k1[2];
      k2[1] = nominal_state[3] + 0.5 * dt * k1[3];
      k2[2] = forward_accel * cos(theta_k2); // Vẫn dùng forward_accel
      k2[3] = forward_accel * sin(theta_k2); // Vẫn dùng forward_accel
      k2[4] = gyro_z_raw;

      // -- Tính k3
      double theta_k3 = nominal_state[4] + 0.5 * dt * k2[4];
      k3[0] = nominal_state[2] + 0.5 * dt * k2[2];
      k3[1] = nominal_state[3] + 0.5 * dt * k2[3];
      k3[2] = forward_accel * cos(theta_k3);
      k3[3] = forward_accel * sin(theta_k3);
      k3[4] = gyro_z_raw;

      // -- Tính k4
      double theta_k4 = nominal_state[4] + dt * k3[4];
      k4[0] = nominal_state[2] + dt * k3[2];
      k4[1] = nominal_state[3] + dt * k3[3];
      k4[2] = forward_accel * cos(theta_k4);
      k4[3] = forward_accel * sin(theta_k4);
      k4[4] = gyro_z_raw;

      // Cập nhật trạng thái
      for (int i = 0; i < NOMINAL_STATE_DIM; i++)
      {
        nominal_state[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
      }
      nominal_state[4] = normalizeAngle(nominal_state[4]);

      // =================== BƯỚC 2: CẬP NHẬT HIỆP PHƯƠNG SAI LỖI (ERROR STATE COVARIANCE) ===================
      double F[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
      double theta_pred = nominal_state[4];
      double cos_th_pred = cos(theta_pred);
      double sin_th_pred = sin(theta_pred);

      ax_raw = accel[0] - error_state[6];
      ay_raw = accel[1] - error_state[7];
      ay_std_body = ax_raw; // Chỉ dùng gia tốc dọc

      for (int i = 0; i < ERROR_STATE_DIM; ++i)
        F[i * ERROR_STATE_DIM + i] = 1.0;

      F[0 * ERROR_STATE_DIM + 2] = dt; // vị trí phụ thuộc trực tiếp vào vận tốc.
      F[1 * ERROR_STATE_DIM + 3] = dt; // vị trí phụ thuộc trực tiếp vào vận tốc.

      // Đạo hàm theo theta (state 4)
      F[2 * ERROR_STATE_DIM + 4] = (-ay_std_body * sin_th_pred) * dt; // vận tốc chịu ảnh hưởng bởi góc quay (ψ).
      F[3 * ERROR_STATE_DIM + 4] = (ay_std_body * cos_th_pred) * dt;  // vận tốc chịu ảnh hưởng bởi góc quay (ψ).

      F[4 * ERROR_STATE_DIM + 5] = -dt; // heading chịu ảnh hưởng bởi bias gyro.

      F[2 * ERROR_STATE_DIM + 6] = -cos_th_pred * dt; // vận tốc chịu ảnh hưởng của bias gia tốc trục x (b_ax).
      F[2 * ERROR_STATE_DIM + 7] = 0;                 // δv_E liên quan đến b_a,y

      // chỗ này có thể sai dấu
      F[3 * ERROR_STATE_DIM + 6] = -sin_th_pred * dt; // vận tốc chịu ảnh hưởng của bias gia tốc trục x (b_ax).
      F[3 * ERROR_STATE_DIM + 7] = 0;                 // δv_N liên quan đến b_a,y

      // double FP[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
      // for (int i = 0; i < ERROR_STATE_DIM; i++)
      // {
      //   for (int j = 0; j < ERROR_STATE_DIM; j++)
      //   {
      //     for (int k = 0; k < ERROR_STATE_DIM; k++)
      //     {
      //       FP[i * ERROR_STATE_DIM + j] += F[i * ERROR_STATE_DIM + k] * P[k * ERROR_STATE_DIM + j];
      //     }
      //   }
      // }

      // Tính FP = F * P, chỉ tính các cột/hàng có phần tử khác 0
      double FP[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        // Chỉ cần tính cho cột 2, 3, 4, 5, 6, 7
        FP[i * ERROR_STATE_DIM + 0] = F[i * ERROR_STATE_DIM + 0] * P[0 * ERROR_STATE_DIM + 0];
        FP[i * ERROR_STATE_DIM + 1] = F[i * ERROR_STATE_DIM + 1] * P[1 * ERROR_STATE_DIM + 1];
        FP[i * ERROR_STATE_DIM + 2] = F[i * ERROR_STATE_DIM + 2] * P[2 * ERROR_STATE_DIM + 2] + F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 2];
        FP[i * ERROR_STATE_DIM + 3] = F[i * ERROR_STATE_DIM + 3] * P[3 * ERROR_STATE_DIM + 3] + F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 3];
        FP[i * ERROR_STATE_DIM + 4] = F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 4] + F[i * ERROR_STATE_DIM + 5] * P[5 * ERROR_STATE_DIM + 4];
        FP[i * ERROR_STATE_DIM + 5] = F[i * ERROR_STATE_DIM + 5] * P[5 * ERROR_STATE_DIM + 5];
        FP[i * ERROR_STATE_DIM + 6] = F[i * ERROR_STATE_DIM + 6] * P[6 * ERROR_STATE_DIM + 6] + F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 6];
        FP[i * ERROR_STATE_DIM + 7] = F[i * ERROR_STATE_DIM + 7] * P[7 * ERROR_STATE_DIM + 7] + F[i * ERROR_STATE_DIM + 4] * P[4 * ERROR_STATE_DIM + 7];
      }

      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = 0; j < ERROR_STATE_DIM; j++)
        {
          double temp = 0;
          for (int k = 0; k < ERROR_STATE_DIM; k++)
          {
            temp += FP[i * ERROR_STATE_DIM + k] * F[j * ERROR_STATE_DIM + k];
          }
          // Thêm nhiễu quá trình Q * dt
          P[i * ERROR_STATE_DIM + j] = temp + Q[i * ERROR_STATE_DIM + j] * dt;
        }
      }

      for (int i = 0; i < ERROR_STATE_DIM; i++)
      {
        for (int j = i + 1; j < ERROR_STATE_DIM; j++)
        {
          P[j * ERROR_STATE_DIM + i] = P[i * ERROR_STATE_DIM + j];
        }
      }

      // for (int i = 0; i < ERROR_STATE_DIM; i++)
      // {
      //   if (P[i * ERROR_STATE_DIM + i] < 0)
      //   {
      //     Serial.print(i);
      //     Serial.println(" Warning: P has negative diagonal elements!");
      //     P[i * ERROR_STATE_DIM + i] = 1e-6; // Đặt lại giá trị nhỏ dương
      //   }
      // }

      // =================== BƯỚC 3 (MỚI): MINI-UPDATE DÙNG IMU YAW ===================
      // prediction_counter++;
      // if (prediction_counter % 10 == 0) // Thực hiện ở tần số 10Hz (100Hz / 10)
      // {
      //   // Ma trận H cho phép đo yaw (1x8)
      //   // H_yaw = [0, 0, 0, 0, 1, 0, 0, 0]

      //   // S_yaw = H_yaw * P * H_yaw^T + R_yaw  (S là số vô hướng 1x1)
      //   double S_yaw = P[4 * ERROR_STATE_DIM + 4] + R_yaw;

      //   // K_yaw = P * H_yaw^T * S_yaw_inv (K là véc-tơ 8x1)
      //   double K_yaw[ERROR_STATE_DIM];
      //   for (int i = 0; i < ERROR_STATE_DIM; i++)
      //   {
      //     K_yaw[i] = P[i * ERROR_STATE_DIM + 4] / S_yaw;
      //   }

      //   // innovation_yaw = measured_yaw - state[4]
      //   // QUAN TRỌNG: Phải chuẩn hóa sai số góc để xử lý việc vượt qua mốc -180/180
      //   double innovation_yaw = normalizeAngle(measured_yaw - nominal_state[4]);

      //   // Cập nhật trạng thái lỗi
      //   for (int i = 0; i < ERROR_STATE_DIM; i++)
      //   {
      //     error_state[i] += K_yaw[i] * innovation_yaw;
      //   }

      //   // Tiêm sai số vào trạng thái danh nghĩa
      //   nominal_state[0] += error_state[0]; // p_E
      //   nominal_state[1] += error_state[1]; // p_N
      //   nominal_state[2] += error_state[2]; // v_E
      //   nominal_state[3] += error_state[3]; // v_N
      //   nominal_state[4] += error_state[4]; // ψ
      //   nominal_state[4] = normalizeAngle(nominal_state[4]);

      //   // Đặt lại trạng thái lỗi
      //   error_state[0] = 0;
      //   error_state[1] = 0;
      //   error_state[2] = 0;
      //   error_state[3] = 0;
      //   error_state[4] = 0;

      //   // Cập nhật hiệp phương sai: P = (I - K * H) * P
      //   double I_KH[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
      //   for (int i = 0; i < ERROR_STATE_DIM; i++)
      //   {
      //     for (int j = 0; j < ERROR_STATE_DIM; j++)
      //     {
      //       // H_yaw chỉ có phần tử thứ 4 là khác 0
      //       double kh_ij = K_yaw[i] * ((j == 4) ? 1.0 : 0.0);
      //       I_KH[i * ERROR_STATE_DIM + j] = ((i == j) ? 1.0 : 0.0) - kh_ij;
      //     }
      //   }

      //   double P_new[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
      //   for (int i = 0; i < ERROR_STATE_DIM; i++)
      //   {
      //     for (int j = 0; j < ERROR_STATE_DIM; j++)
      //     {
      //       for (int k = 0; k < ERROR_STATE_DIM; k++)
      //       {
      //         P_new[i * ERROR_STATE_DIM + j] += I_KH[i * ERROR_STATE_DIM + k] * P[k * ERROR_STATE_DIM + j];
      //       }
      //     }
      //   }
      //   memcpy(P, P_new, sizeof(P_new));

      //   for (int i = 0; i < ERROR_STATE_DIM; i++)
      //   {
      //     for (int j = i + 1; j < ERROR_STATE_DIM; j++)
      //     {
      //       P[j * ERROR_STATE_DIM + i] = P[i * ERROR_STATE_DIM + j];
      //     }
      //   }

      //   for (int i = 0; i < ERROR_STATE_DIM; i++)
      //   {
      //     if (P[i * ERROR_STATE_DIM + i] < 0)
      //     {
      //       Serial.print(i);
      //       Serial.println(" Warning: P has negative diagonal elements!");
      //       P[i * ERROR_STATE_DIM + i] = 1e-6; // Đặt lại giá trị nhỏ dương
      //     }
      //   }
      // }
      // Sau mỗi bước, in ra state
      // Serial.print("x=");
      // Serial.print(nominal_state[0], 3);
      // Serial.print("  y=");
      // Serial.print(nominal_state[1], 3);
      // Serial.print("  vx=");
      // Serial.print(nominal_state[2], 3);
      // Serial.print("  vy=");
      // Serial.print(nominal_state[3], 3);
      // Serial.print("  yaw=");
      // Serial.print(nominal_state[4] * 180.0 / M_PI, 2);
      // Serial.print("  errbax=");
      // Serial.print(error_state[6], 3);
      // Serial.print("  errbay=");
      // Serial.print(error_state[7], 3);
      // Serial.print("  yaw=");
      // Serial.println(nominal_state[4] * 180.0 / M_PI, 2);
      xSemaphoreGive(stateMutex);
    }
    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
}

// Tính determinant của ma trận n x n (dùng đệ quy, phù hợp n <= 5)
double computeDeterminant(double *M, int n)
{
  if (n == 1)
    return M[0];
  if (n == 2)
    return M[0] * M[3] - M[1] * M[2];

  double det = 0.0;
  double sub[25]; // đủ chứa ma trận con tối đa 5x5
  for (int x = 0; x < n; x++)
  {
    int subi = 0;
    for (int i = 1; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (j == x)
          continue;
        sub[subi++] = M[i * n + j];
      }
    }
    det += ((x % 2 == 0 ? 1 : -1) * M[0 * n + x] * computeDeterminant(sub, n - 1));
  }
  return det;
}

// Hàm cho riêng 5x5
double computeDeterminant5x5(const double S[5 * 5])
{
  return computeDeterminant((double *)S, 5);
}

void updateTask(void *pvParameters)
{
  const TickType_t xFrequency = 100 / portTICK_PERIOD_MS; // 10Hz
  TickType_t xLastWakeTime = xTaskGetTickCount();
  uint8_t rtksta = 0;
  int32_t hAcc = 0;
  int32_t vAcc = 0;
  int32_t headingAcc = 0;
  PointECEF tempECEF;
  PointECEF tempECEF_2;
  while (1)
  {
    bool PVT = false;

    double z[5] = {0, 0, 0, 0, 0}; // Vẫn là [current_east, current_north, current_east_velocity, current_north_velocity, current_yaw]
    tempECEF.x = 0;
    tempECEF.y = 0;
    tempECEF.z = 0;
    // if (myGNSS.getNAVHPPOSECEF())
    if (myGNSS.getPVT())
    {
      // rtksta = myGNSS.getCarrierSolutionType();
      PointLLH tempLLH;

      // Lấy tọa độ và chuyển sang hệ ECEF
      // GetRTK(tempLLH);
      // Serial.print(tempLLH.latitude, 9);
      // Serial.print(", ");
      // Serial.print(tempLLH.longitude, 9);
      // Serial.print(", ");
      // llhToEcef(tempLLH, tempECEF_2);
      // First, let's collect the position data
      int32_t ECEFX = myGNSS.getHighResECEFX();
      int8_t ECEFXHp = myGNSS.getHighResECEFXHp();
      int32_t ECEFY = myGNSS.getHighResECEFY();
      int8_t ECEFYHp = myGNSS.getHighResECEFYHp();
      int32_t ECEFZ = myGNSS.getHighResECEFZ();
      int8_t ECEFZHp = myGNSS.getHighResECEFZHp();
      uint32_t accuracy = myGNSS.getPositionAccuracy();

      // Defines storage for the ECEF coordinates as double
      double d_ECEFX;
      double d_ECEFY;
      double d_ECEFZ;

      // Assemble the high precision coordinates
      d_ECEFX = ((double)ECEFX) / 100.0;      // Convert from cm to m
      d_ECEFX += ((double)ECEFXHp) / 10000.0; // Now add the high resolution component ( mm * 10^-1 = m * 10^-4 )
      d_ECEFY = ((double)ECEFY) / 100.0;      // Convert from cm to m
      d_ECEFY += ((double)ECEFYHp) / 10000.0; // Now add the high resolution component ( mm * 10^-1 = m * 10^-4 )
      d_ECEFZ = ((double)ECEFZ) / 100.0;      // Convert from cm to m
      d_ECEFZ += ((double)ECEFZHp) / 10000.0; // Now add the high resolution component ( mm * 10^-1 = m * 10^-4 )

      tempECEF.x = d_ECEFX;
      tempECEF.y = d_ECEFY;
      tempECEF.z = d_ECEFZ;

      // BƯỚC 2: Chuyển đổi tọa độ ECEF vừa nhận được sang ENU
      double current_east, current_north, current_up;
      ecefToEnu(tempECEF, current_east, current_north, current_up);

      double vel_east = ((double)myGNSS.getNedEastVel()) / 1000.0;   // Chuyển từ mm/s sang m/s
      double vel_north = ((double)myGNSS.getNedNorthVel()) / 1000.0; // Chuyển từ mm/s sang m/s
      double heading = ((double)myGNSS.getHeading()) / 100000.0;     // Chuyển từ 1e-5 độ sang độ
      // heading = 0.0;
      // current_east = 1.0;
      // current_north = 2.0;
      // vel_east = 0.0;
      // vel_north = 0.0;
      // heading = 1.0;

      z[0] = current_east;
      z[1] = current_north;
      z[2] = vel_east;
      z[3] = vel_north;
      z[4] = heading * M_PI / 180.0; // Chuyển sang radian
      z[4] = normalizeAngle(z[4]);
      // z[4] = nominal_state[4]; // Chỉ cập nhật vị trí, không cập nhật yaw

      //*********DEBUG***********//
      // Serial.print(current_east_2,4);
      // Serial.print(", ");
      // Serial.print(current_north_2,4);
      // Serial.print(" | ");
      // Serial.print(current_east, 4);
      // Serial.print(", ");
      // Serial.print(current_north, 4);
      // Serial.print(", ");
      // Serial.print(sqrt((tempECEF.x - enu_origin_ecef.x) * (tempECEF.x - enu_origin_ecef.x) + (tempECEF.y - enu_origin_ecef.y) * (tempECEF.y - enu_origin_ecef.y)));
      // Serial.print(", ");
      //*********DEBUG***********//

      rtksta = myGNSS.getCarrierSolutionType();
      hAcc = myGNSS.getHorizontalAccEst();    // mm
      vAcc = myGNSS.getVerticalAccEst();      // mm
      headingAcc = myGNSS.getHeadingAccEst(); // 1e-5 độ
      if (myGNSS.getGroundSpeed() < 300)      // Nếu tốc độ mặt đất < 0.3m/s thì không tin tưởng vào độ chính xác hướng
        headingAcc = 3500000;                 // Đặt độ chính xác hướng rất lớn (180 độ)
      // headingAcc = 10;
      PVT = true; // Đánh dấu đã có dữ liệu hợp lệ
    }

    // if (myGNSS.getPVT() && (myGNSS.getInvalidLlh() == false))
    // {
    //   rtksta = myGNSS.getCarrierSolutionType();
    //   hAcc = myGNSS.getHorizontalAccEst(); // mm
    //   PVT = true;
    // }

    if (PVT)
    {
      if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
      {
        Serial.print(z[0]);
        Serial.print(",");
        Serial.print(z[1]);
        Serial.print(",");
        Serial.print(z[2]);
        Serial.print(",");
        Serial.print(z[3]);
        Serial.print(",");
        Serial.print(z[4] * 180.0 / M_PI);
        Serial.print(",   ");

        Serial.print(nominal_state[0], 2);
        Serial.print(",");
        Serial.print(nominal_state[1], 2);
        Serial.print(",");
        Serial.print(nominal_state[2], 2);
        Serial.print(",");
        Serial.print(nominal_state[3], 2);
        Serial.print(",");
        Serial.print((nominal_state[4] * 180.0 / M_PI), 2);
        Serial.print(",");
        Serial.print(gocIMU, 2);
        Serial.print(",   ");

        // =================== BƯỚC 3: TÍNH TOÁN HIỆU CHỈNH LỖI ===================

        // Ma trận H cho trạng thái lỗi vẫn như cũ
        static double H[5 * ERROR_STATE_DIM] = {0};
        H[0 * ERROR_STATE_DIM + 0] = 1.0; // d(z_e)/d(δp_e)
        H[1 * ERROR_STATE_DIM + 1] = 1.0; // d(z_n)/d(δp_n)
        H[2 * ERROR_STATE_DIM + 2] = 1.0; // d(z_ve)/d(δv_e)
        H[3 * ERROR_STATE_DIM + 3] = 1.0; // d(z_vn)/d(δv_n)
        H[4 * ERROR_STATE_DIM + 4] = 1.0; // d(z_ψ)/d(δψ)

        // =================== BẮT ĐẦU ADAPTIVE R ===================
        double sigma_p = (double)hAcc / 1000.0;                           // mm -> m
        double sigma_v = (double)vAcc / 1000.0;                           // mm/s -> m/s
        double sigma_yaw = ((double)headingAcc / 10000.0) * M_PI / 180.0; // độ -> radian

        // Serial.print("sigma_p:");
        // Serial.print(sigma_p, 4);
        // Serial.print(", sigma_v:");
        // Serial.print(sigma_v, 4);
        // Serial.print(", sigma_yaw:");
        // Serial.print(sigma_yaw * 180.0 / M_PI, 4);

        static double R[5 * 5] = {0};
        R[0 * 5 + 0] = sigma_p * sigma_p;
        R[1 * 5 + 1] = sigma_p * sigma_p;
        R[2 * 5 + 2] = sigma_v * sigma_v;
        R[3 * 5 + 3] = sigma_v * sigma_v;
        R[4 * 5 + 4] = sigma_yaw * sigma_yaw;
        // Serial.print(R_local[0], 5);
        // Serial.print(", ");
        // =================== KẾT THÚC ADAPTIVE R ===================

        // =================== TÍNH TOÁN TỐI ƯU ===================
        // --- BƯỚC 1: TÍNH S = H*P*H^T + R ---
        // Tính S = H * P * H^T + R

        // for (int i = 0; i < ERROR_STATE_DIM; i++)
        // {
        //   if (P[i * ERROR_STATE_DIM + i] < 1e-6)
        //   {
        //     Serial.print("Warning: P[");
        //     Serial.print(i);
        //     Serial.println("][");
        //     Serial.print(i);
        //     Serial.println("] too small!");
        //     P[i * ERROR_STATE_DIM + i] = 1e-6;
        //   }
        // }

        // Tính H_P = H * P
        double H_P[5 * ERROR_STATE_DIM] = {0};
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

        // Tính S = H_P * H^T + R
        static double S[5 * 5] = {0};
        for (int i = 0; i < 5; i++)
        {
          for (int j = 0; j < 5; j++)
          {
            for (int k = 0; k < ERROR_STATE_DIM; k++)
            {
              S[i * 5 + j] += H_P[i * ERROR_STATE_DIM + k] * H[j * ERROR_STATE_DIM + k]; // Sửa H[k * ERROR_STATE_DIM + j] thành H[j * ERROR_STATE_DIM + k]
            }
            if (i == j)
              S[i * 5 + j] += R[i * 5 + j];
          }
        }

        // Kiểm tra S khả nghịch
        double det = computeDeterminant5x5(S); // Hàm tính định thức
        if (fabs(det) < 1e-6)
        {
          Serial.println("Warning: Cannot invert S!");
          xSemaphoreGive(stateMutex);
          vTaskDelayUntil(&xLastWakeTime, xFrequency);
          continue;
        }

        // Tính nghịch đảo S (5x5)
        static double S_inv[5 * 5];
        if (!inverse_5x5(S, S_inv))
        {
          Serial.println("Warning: Matrix inversion failed!");
          xSemaphoreGive(stateMutex);
          vTaskDelayUntil(&xLastWakeTime, xFrequency);
          continue;
        }

        // Tính Kalman gain K = P * H^T * S^-1
        double P_HT[ERROR_STATE_DIM * 5] = {0};
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
        double K[ERROR_STATE_DIM * 5] = {0};
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

        // --- BƯỚC 3: CẬP NHẬT TRẠNG THÁI ---
        // Tính innovation
        double innovation[5] = {
            z[0] - nominal_state[0],
            z[1] - nominal_state[1],
            z[2] - nominal_state[2],
            z[3] - nominal_state[3],
            normalizeAngle(z[4] - nominal_state[4])};

        // Cập nhật trạng thái lỗi
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          for (int j = 0; j < 5; j++)
          {
            error_state[i] += K[i * 5 + j] * innovation[j];
          }
        }

        // Serial.print("|error:");
        // Serial.print(error_state[0], 2);
        // Serial.print(",");
        // Serial.print(error_state[1], 2);
        // Serial.print(",");
        // Serial.print(error_state[2], 2);
        // Serial.print(",");
        // Serial.print(error_state[3], 2);
        // Serial.print(",");
        // Serial.print(error_state[4] * 180.0 / M_PI, 2);
        // Serial.print(",");
        // Serial.print(error_state[5] * 180.0 / M_PI, 2);
        // Serial.print(",");
        // Serial.print(error_state[6], 2);
        // Serial.print(",");
        // Serial.print(error_state[7], 2);
        // Serial.print(" | ");

        // Tiêm lỗi vào trạng thái danh nghĩa
        nominal_state[0] += error_state[0]; // p_E
        nominal_state[1] += error_state[1]; // p_N
        nominal_state[2] += error_state[2]; // v_E
        nominal_state[3] += error_state[3]; // v_N
        nominal_state[4] += error_state[4]; // ψ
        nominal_state[4] = normalizeAngle(nominal_state[4]);

        // Đặt lại trạng thái lỗi
        error_state[0] = 0; // δp_E
        error_state[1] = 0; // δp_N
        error_state[2] = 0; // δv_E
        error_state[3] = 0; // δv_N
        error_state[4] = 0; // δψ
        // Giữ nguyên b_g,z, b_a,x, b_a,y

        // Cập nhật ma trận hiệp phương sai P = (I - K*H) * P
        static double I_minus_KH[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          for (int j = 0; j < ERROR_STATE_DIM; j++)
          {
            I_minus_KH[i * ERROR_STATE_DIM + j] = (i == j) ? 1.0 : 0.0;
            for (int k = 0; k < 5; k++)
            {
              I_minus_KH[i * ERROR_STATE_DIM + j] -= K[i * 5 + k] * H[k * ERROR_STATE_DIM + j];
            }
          }
        }
        static double P_new[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
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
        for (int i = 0; i < ERROR_STATE_DIM; i++)
        {
          for (int j = i + 1; j < ERROR_STATE_DIM; j++)
          {
            P[j * ERROR_STATE_DIM + i] = P[i * ERROR_STATE_DIM + j];
          }
        }

        // =================== BƯỚC 4: MINI-UPDATE DÙNG IMU YAW (RELATIVE ANGLE) ===================
        // Mục đích: Dùng độ mượt của IMU Yaw để hiệu chỉnh độ trôi (b_gz) và làm mượt dpsi
        if (fabs(C0_yaw) > 1e-4) // Chỉ thực hiện nếu C0_yaw đã được tính toán ở setup
        {
          // 1. Góc IMU đã được chuyển sang radian (gocIMU đang là độ)
          double measured_imu_yaw = gocIMU * M_PI / 180.0;

          // 2. Áp dụng góc bù trừ C0_yaw
          // Đây là góc IMU đã được chuyển về hệ Navigation (North-referenced)
          double corrected_imu_yaw = normalizeAngle(measured_imu_yaw - C0_yaw);

          // 3. Tính toán Innovation (Sai số đo lường)
          // Sai số = corrected_imu_yaw - nominal_state[4] (yaw hiện tại)
          double innovation_imu_yaw = normalizeAngle(corrected_imu_yaw - nominal_state[4]);

          // Ma trận H cho phép đo yaw (1x8)
          // H_yaw = [0, 0, 0, 0, 1, 0, 0, 0]
          // Ma trận R cho phép đo yaw (1x1) là R_yaw (đã được tính ở setup, R_yaw nhỏ)

          // S_yaw = H_yaw * P * H_yaw^T + R_yaw  (S là số vô hướng 1x1)
          double S_yaw = P[4 * ERROR_STATE_DIM + 4] + R_yaw;

          // K_yaw = P * H_yaw^T * S_yaw_inv (K là véc-tơ 8x1)
          double K_yaw[ERROR_STATE_DIM];
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
          double I_KH[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
          for (int i = 0; i < ERROR_STATE_DIM; i++)
          {
            for (int j = 0; j < ERROR_STATE_DIM; j++)
            {
              // H_yaw chỉ có phần tử thứ 4 là khác 0 (trạng thái dpsi)
              double kh_ij = K_yaw[i] * ((j == 4) ? 1.0 : 0.0);
              I_KH[i * ERROR_STATE_DIM + j] = ((i == j) ? 1.0 : 0.0) - kh_ij;
            }
          }

          double P_new_imu[ERROR_STATE_DIM * ERROR_STATE_DIM] = {0};
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
        } // Kết thúc Update IMU

        // Trạng thái lỗi đã được "tiêu thụ" và được reset về 0 cho vòng lặp tiếp theo.
        // Hiệp phương sai P đã được cập nhật để phản ánh sự không chắc chắn mới.
        Serial.print(rtksta);
        Serial.print(",");
        Serial.print(hAcc);
        Serial.print(",");

        Serial.print(nominal_state[0], 3);
        Serial.print(",");
        Serial.print(nominal_state[1], 3);
        Serial.print(",");
        Serial.print(nominal_state[2], 3);
        Serial.print(",");
        Serial.print(nominal_state[3], 3);
        Serial.print(",");
        Serial.print(nominal_state[4], 3);
        Serial.print(",");
        Serial.print(nominal_state[4] * 180 / M_PI, 3);
        Serial.println();
        xSemaphoreGive(stateMutex);
      }
    }
  }
  vTaskDelayUntil(&xLastWakeTime, xFrequency);
}

static void SensorUartSend(uint8_t *p_data, uint32_t uiSize);
static void CopeSensorData(uint32_t uiReg, uint32_t uiRegNum);
static void Delayms(uint16_t ucMs);

static void SensorUartSend(uint8_t *p_data, uint32_t uiSize)
{
  Serial1.write(p_data, uiSize);
  Serial1.flush();
}

static void CopeSensorData(uint32_t uiReg, uint32_t uiRegNum)
{
  int i;
  for (i = 0; i < uiRegNum; i++)
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
void networkTask(void *parameter)
{
  //   unsigned long lastSendTime = 0;
  //   ESKF_STATE_DATA datasend;
  //   while (1)
  //   {
  //     if (SerialBT.connected())
  //     {
  //       // Nhận Tham số Tinh chỉnh (TUNING_CONFIG)
  //       if (SerialBT.available() >= sizeof(TUNING_CONFIG))
  //       {
  //         TUNING_CONFIG newConfig;
  //         size_t bytesRead = SerialBT.readBytes((uint8_t *)&newConfig, sizeof(TUNING_CONFIG));
  //         if (bytesRead == sizeof(TUNING_CONFIG))
  //         {
  //           if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
  //           {
  //             currentConfig = newConfig;
  //             Serial.print("Received Q: ");
  //             for (int i = 0; i < 3; i++)
  //               Serial.printf("%.2f ", currentConfig.Q_noise_params[i]);
  //             Serial.print("| R: ");
  //             for (int i = 0; i < 3; i++)
  //               Serial.printf("%.2f ", currentConfig.R_noise_params[i]);
  //             Serial.println(". Updated ESKF.");
  //             xSemaphoreGive(stateMutex);
  //           }
  //         }
  //       }

  //       // Gửi Dữ liệu Trạng thái (ESKF_STATE_DATA)
  //       if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
  //       {
  //         datasend.timestamp = millis();
  //         datasend.state_x[0] = 1.0f + (millis() / 1000.0);
  //         datasend.state_x[1] = 2.0f + (millis() / 1000.0);
  //         datasend.state_x[2] = 3.0f + (millis() / 1000.0);
  //         datasend.state_x[3] = 4.0f + (millis() / 1000.0);
  //         datasend.state_x[4] = 5.0f + (millis() / 1000.0);
  //         for (int i = 0; i < ERROR_STATE_DIM; i++)
  //         {
  //           datasend.P_diag[i] = 0.2f + (i * 0.1f * (millis() / 1000.0));
  //         }
  //         datasend.raw_data[0] = 9.81f + (millis() / 1000.0);
  //         datasend.raw_data[1] = 0.1f + (millis() / 1000.0);
  //         datasend.raw_data[2] = 0.2f + (millis() / 1000.0);
  //         datasend.raw_data[3] = 0.2f + (millis() / 1000.0);
  //         datasend.raw_data[4] = 0.2f + (millis() / 1000.0);
  //         datasend.raw_data[5] = 0.2f + (millis() / 1000.0);
  //         datasend.raw_data[6] = 0.2f + (millis() / 1000.0);
  //         datasend.raw_data[7] = 0.2f + (millis() / 1000.0);
  //         datasend.raw_data[8] = 0.2f + (millis() / 1000.0);
  //         datasend.gnss_accuracy[0] = 0.2f + (millis() / 1000.0);
  //         datasend.gnss_accuracy[1] = 0.2f + (millis() / 1000.0);
  //         datasend.gnss_accuracy[2] = 0.2f + (millis() / 1000.0);

  //         size_t bytesSent = SerialBT.write((uint8_t *)&datasend, sizeof(ESKF_STATE_DATA));
  //         if (bytesSent == sizeof(ESKF_STATE_DATA))
  //         {
  //           Serial.println("Data sent successfully.");
  //         }
  //         else
  //         {
  //           Serial.println("Failed to send data.");
  //         }
  //         xSemaphoreGive(stateMutex);
  //       }
  //     }
  //     else
  //     {
  //       Serial.println("Waiting for Bluetooth connection...");
  //     }

  //     vTaskDelay(pdMS_TO_TICKS(20)); // Giảm xuống 50Hz để tránh nghẽn Bluetooth
  //   }
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
        CalculateEnuMatrix(initialLLH);

        // Lấy tọa độ gốc ECEF trực tiếp từ GNSS
        enu_origin_ecef.x = ((double)myGNSS.getHighResECEFX() / 100.0) + ((double)myGNSS.getHighResECEFXHp() / 10000.0);
        enu_origin_ecef.y = ((double)myGNSS.getHighResECEFY() / 100.0) + ((double)myGNSS.getHighResECEFYHp() / 10000.0);
        enu_origin_ecef.z = ((double)myGNSS.getHighResECEFZ() / 100.0) + ((double)myGNSS.getHighResECEFZHp() / 10000.0);

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

    double sigma_p = (double)hAcc_p / 1000.0;                    // mm -> m
    double sigma_v = vAcc_p > 0 ? (double)vAcc_p / 1000.0 : 0.1; // mm/s -> m/s, mặc định 0.1 m/s nếu không có
    double sigma_yaw;                                            // độ -> radian, mặc định 5°
    double sigma_bg_z = 0.1;                                     // rad/s (bias con quay)
    double sigma_ba = 0.3;
    bool heading_valid = false;
    while (!heading_valid)
    {
      if (myGNSS.getGroundSpeed() > V_MIN_HEADING * 1000)
      {
        heading_valid = true;
        headingAcc_p = myGNSS.getHeadingAccEst();
        sigma_yaw = headingAcc_p > 0 ? (double)headingAcc_p / 10000.0 * M_PI / 180.0 : 0.0873; // độ -> radian, mặc định 5°
        // sigma_yaw = 5.0 * M_PI / 180.0; // Đặt độ chính xác hướng rất lớn (5 độ)
        nominal_state[4] = myGNSS.getHeading() / 100000.0 * M_PI / 180.0;
        nominal_state[4] = normalizeAngle(nominal_state[4]);
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
    Serial.print(sigma_yaw * 180.0 / M_PI, 4);

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
    double initial_imu_yaw_deg = sReg[Yaw] / 32768.0f * 180.0f;
    C0_yaw = initial_imu_yaw_deg * M_PI / 180.0 - nominal_state[4];
    s_cDataUpdate &= ~ANGLE_UPDATE;
    Serial.print("Initial IMU yaw offset (C0_yaw) calculated: ");
    Serial.println(C0_yaw * 180.0 / M_PI);
  }
  else
  {
    C0_yaw = 0.0;
    Serial.println("Warning: Could not read initial IMU heading.");
  }

  // --- 5.4: Khởi tạo ma trận nhiễu quá trình Q ---
  double vel_noise_std_dev = 1.5;
  double accel_bias_noise_std_dev = 0.86999;
  double gyro_bias_noise_std_dev = 0.8994;

  Q[2 * ERROR_STATE_DIM + 2] = vel_noise_std_dev * vel_noise_std_dev;
  Q[3 * ERROR_STATE_DIM + 3] = vel_noise_std_dev * vel_noise_std_dev;
  Q[4 * ERROR_STATE_DIM + 4] = (100.0 * M_PI / 180.0) * (100.0 * M_PI / 180.0); // Nhiễu đo của Gyro
  Q[5 * ERROR_STATE_DIM + 5] = gyro_bias_noise_std_dev * gyro_bias_noise_std_dev;
  Q[6 * ERROR_STATE_DIM + 6] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;
  Q[7 * ERROR_STATE_DIM + 7] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;

  // Khởi tạo ma trận nhiễu đo lường R cho IMU Yaw ---
  double yaw_measurement_std_dev_rad = 0.5 * M_PI / 180.0; // Sai số đo của IMU yaw là 1.1 độ
  R_yaw = yaw_measurement_std_dev_rad * yaw_measurement_std_dev_rad;

  // KHỐI 6: TẠO CÁC TÁC VỤ (TASKS)
  Serial.println("Creating FreeRTOS tasks...");
  xTaskCreatePinnedToCore(predictionTask, "Prediction", 8192, NULL, 5, NULL, 0); // Core 0
  xTaskCreatePinnedToCore(updateTask, "Update", 12288, NULL, 5, NULL, 1);        // Core 1
  Serial.println("System initialized successfully. ES-EKF is running.");
  delay(2000);
}
  // if (s_cDataUpdate)
  // {
  //   for (i = 0; i < 3; i++)
  //   {
  //     fAcc[i] = sReg[AX + i] / 32768.0f * 16.0f;
  //     fGyro[i] = sReg[GX + i] / 32768.0f * 2000.0f;
  //     fAngle[i] = sReg[Roll + i] / 32768.0f * 180.0f;
  //     double initial_imu_yaw_deg = sReg[Yaw] / 32768.0f * 180.0f;
  //   }
  //   if (s_cDataUpdate & ANGLE_UPDATE)
  //   {
  //     double initial_imu_yaw_deg = fAngle[2];
  //     C0_yaw = initial_imu_yaw_deg * M_PI / 180.0 - nominal_state[4];
  //     // nominal_state[4] = fAngle[2];
  //     s_cDataUpdate &= ~ANGLE_UPDATE;
  //   }
  //   s_cDataUpdate = 0;
  // }
  // else
  // {
  //   nominal_state[4] = 0; // Mặc định 0 nếu không đọc được
  //   Serial.println("Warning: Initial heading not available, set to 0.");
  // }
  // // nominal_state[4] = 0;
  // //  if (initialHeading.isDataValid != 0)
  // //  {                                                       // Kiểm tra dữ liệu hợp lệ
  // //    nominal_state[4] = initialHeading.yaw * M_PI / 180.0; // Chuyển sang radian
  // //  }

  // Serial.print("Initial Heading (rad): ");
  // Serial.println(nominal_state[4]);

  // if (stateMutex == NULL)
  // {
  //   Serial.println("Error: Failed to create mutex!");
  //   while (1)
  //     ; // Treo hệ thống
  // }
  // double yaw_measurement_std_dev_rad = 1.1 * M_PI / 180.0;
  // R_yaw = yaw_measurement_std_dev_rad * yaw_measurement_std_dev_rad;

  // // =================== BỘ THAM SỐ ĐỀ XUẤT ===================

  // // Khởi tạo Q - Mức độ nhiễu của quá trình THEO THỜI GIAN
  // // double vel_noise_std_dev = 0.05;             // Giả định vận tốc có thể nhiễu ngẫu nhiên 5 cm/s
  // // double accel_bias_noise_std_dev = 0.026999; // Giả định bias gia tốc RẤT ỔN ĐỊNH, ít thay đổi
  // // double gyro_bias_noise_std_dev = 0.01994;    // Giả định bias gyro CỰC KỲ ỔN ĐỊNH
  // double vel_noise_std_dev = 1.5;            // Giả định vận tốc có thể nhiễu ngẫu nhiên 5 cm/s
  // double accel_bias_noise_std_dev = 0.86999; // Giả định bias gia tốc RẤT ỔN ĐỊNH, ít thay đổi
  // double gyro_bias_noise_std_dev = 0.8994;   // Giả định bias gyro CỰC KỲ ỔN ĐỊNH

  // Q[0 * ERROR_STATE_DIM + 0] = 0.0; // δp_E
  // Q[1 * ERROR_STATE_DIM + 1] = 0.0; // δp_N
  // Q[2 * ERROR_STATE_DIM + 2] = vel_noise_std_dev * vel_noise_std_dev;
  // Q[3 * ERROR_STATE_DIM + 3] = vel_noise_std_dev * vel_noise_std_dev;
  // Q[4 * ERROR_STATE_DIM + 4] = (100.0 * M_PI / 180.0) * (100.0 * M_PI / 180.0); // Nhiễu đo của Gyro (từ datasheet)
  // Q[5 * ERROR_STATE_DIM + 5] = gyro_bias_noise_std_dev * gyro_bias_noise_std_dev;
  // Q[6 * ERROR_STATE_DIM + 6] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;
  // Q[7 * ERROR_STATE_DIM + 7] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;

  // // Tạo tác vụ cho hai core

  // xTaskCreatePinnedToCore(predictionTask, "Prediction", 4096, NULL, 5, NULL, 0); // Core 0
  // xTaskCreatePinnedToCore(updateTask, "Update", 4096, NULL, 5, NULL, 1);         // Core 1
  // // xTaskCreatePinnedToCore(networkTask, "Network", 4096, NULL, 3, NULL, 1);

  // Serial.println("System initialized with RTK Fixed position and WTGAHRS3 heading.");
  // delay(2000);
//}

void loop() {}