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

// --- CẤU HÌNH BỘ LỌC ---
const int STATE_DIM = 5; // Kích thước của cả Nominal State và Error State

// --- BIẾN TOÀN CỤC ---
SFE_UBLOX_GNSS myGNSS;

esp_now_peer_info_t peerInfo;
uint8_t broadcastAddress[] = {0x24, 0xD7, 0xEB, 0x13, 0x74, 0xAC};

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
double nominal_state[STATE_DIM] = {0}; // [p_e, p_n, v_e, v_n, theta, b_ax, b_ay, b_gz]

// Ma trận hiệp phương sai của TRẠNG THÁI LỖI (Error State Covariance)
double P[STATE_DIM * STATE_DIM] = {0};

// Ma trận nhiễu quá trình cho TRẠNG THÁI LỖI
double Q[STATE_DIM * STATE_DIM] = {0};

double R_yaw;

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
        accel[0] = 0.0;
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

      double ax_raw = accel[0] - nominal_state[5];
      double ay_raw = accel[1] - nominal_state[6];
      double gyro_z_raw = gyro_z - nominal_state[7];

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
      // Serial.print("ax_nav_east=");
      // Serial.print(ax_nav_east, 3);
      // Serial.print("  ay_nav_north=");
      // Serial.print(ay_nav_north, 3);
      double k1[STATE_DIM], k2[STATE_DIM], k3[STATE_DIM], k4[STATE_DIM];

      // -- Tính k1
      k1[0] = nominal_state[2];
      k1[1] = nominal_state[3];
      k1[2] = ax_nav_east;
      k1[3] = ay_nav_north;
      k1[4] = gyro_z_raw;
      k1[5] = 0;
      k1[6] = 0;
      k1[7] = 0;

      // --Tính k2
      double theta_k2 = nominal_state[4] + 0.5 * dt * k1[4];
      k2[0] = nominal_state[2] + 0.5 * dt * k1[2];
      k2[1] = nominal_state[3] + 0.5 * dt * k1[3];
      k2[2] = forward_accel * cos(theta_k2); // Vẫn dùng forward_accel
      k2[3] = forward_accel * sin(theta_k2); // Vẫn dùng forward_accel
      k2[4] = gyro_z_raw;
      for (int i = 5; i < STATE_DIM; ++i)
        k2[i] = 0;

      // -- Tính k3
      double theta_k3 = nominal_state[4] + 0.5 * dt * k2[4];
      k3[0] = nominal_state[2] + 0.5 * dt * k2[2];
      k3[1] = nominal_state[3] + 0.5 * dt * k2[3];
      k3[2] = forward_accel * cos(theta_k3);
      k3[3] = forward_accel * sin(theta_k3);
      k3[4] = gyro_z_raw;
      for (int i = 5; i < STATE_DIM; ++i)
        k3[i] = 0;

      // -- Tính k4
      double theta_k4 = nominal_state[4] + dt * k3[4];
      k4[0] = nominal_state[2] + dt * k3[2];
      k4[1] = nominal_state[3] + dt * k3[3];
      k4[2] = forward_accel * cos(theta_k4);
      k4[3] = forward_accel * sin(theta_k4);
      k4[4] = gyro_z_raw;
      for (int i = 5; i < STATE_DIM; ++i)
        k4[i] = 0;

      // Cập nhật trạng thái
      for (int i = 0; i < STATE_DIM; i++)
      {
        nominal_state[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
      }
      nominal_state[4] = normalizeAngle(nominal_state[4]);

      // =================== BƯỚC 2: CẬP NHẬT HIỆP PHƯƠNG SAI LỖI (ERROR STATE COVARIANCE) ===================
      double F[STATE_DIM * STATE_DIM] = {0};
      double theta_pred = nominal_state[4];
      double cos_th_pred = cos(theta_pred);
      double sin_th_pred = sin(theta_pred);

      ax_raw = accel[0] - nominal_state[5];
      ay_std_body = ax_raw; // Chỉ dùng gia tốc dọc

      for (int i = 0; i < STATE_DIM; ++i)
        F[i * STATE_DIM + i] = 1.0;

      F[0 * STATE_DIM + 2] = dt;
      F[1 * STATE_DIM + 3] = dt;

      // Đạo hàm theo theta (state 4)
      F[2 * STATE_DIM + 4] = (-ay_std_body * sin_th_pred) * dt;
      F[3 * STATE_DIM + 4] = (ay_std_body * cos_th_pred) * dt;

      // Đạo hàm theo bias_ax (state 5)
      F[2 * STATE_DIM + 5] = -(-sin_th_pred) * dt; // = sin_th_pred * dt
      F[3 * STATE_DIM + 5] = (-cos_th_pred) * dt;  // = -cos_th_pred * dt

      // Đạo hàm theo bias_ay (state 6)
      F[2 * STATE_DIM + 6] = -sin_th_pred * dt;
      F[3 * STATE_DIM + 6] = cos_th_pred * dt;

      F[4 * STATE_DIM + 7] = -dt;

      double FP[STATE_DIM * STATE_DIM] = {0};
      for (int i = 0; i < STATE_DIM; i++)
      {
        for (int j = 0; j < STATE_DIM; j++)
        {
          for (int k = 0; k < STATE_DIM; k++)
          {
            FP[i * STATE_DIM + j] += F[i * STATE_DIM + k] * P[k * STATE_DIM + j];
          }
        }
      }

      for (int i = 0; i < STATE_DIM; i++)
      {
        for (int j = 0; j < STATE_DIM; j++)
        {
          double temp = 0;
          for (int k = 0; k < STATE_DIM; k++)
          {
            temp += FP[i * STATE_DIM + k] * F[j * STATE_DIM + k];
          }
          // Thêm nhiễu quá trình Q * dt
          P[i * STATE_DIM + j] = temp + Q[i * STATE_DIM + j] * dt;
        }
      }
      // =================== BƯỚC 3 (MỚI): MINI-UPDATE DÙNG IMU YAW ===================
      prediction_counter++;
      if (prediction_counter % 10 == 0) // Thực hiện ở tần số 10Hz (100Hz / 10)
      {
        // Ma trận H cho phép đo yaw (1x8)
        // H_yaw = [0, 0, 0, 0, 1, 0, 0, 0]

        // S_yaw = H_yaw * P * H_yaw^T + R_yaw  (S là số vô hướng 1x1)
        double S_yaw = P[4 * STATE_DIM + 4] + R_yaw;

        // K_yaw = P * H_yaw^T * S_yaw_inv (K là véc-tơ 8x1)
        double K_yaw[STATE_DIM];
        for (int i = 0; i < STATE_DIM; i++)
        {
          K_yaw[i] = P[i * STATE_DIM + 4] / S_yaw;
        }

        // innovation_yaw = measured_yaw - state[4]
        // QUAN TRỌNG: Phải chuẩn hóa sai số góc để xử lý việc vượt qua mốc -180/180
        double innovation_yaw = normalizeAngle(measured_yaw - nominal_state[4]);

        // Cập nhật trạng thái: state = state + K * innovation
        for (int i = 0; i < STATE_DIM; i++)
        {
          nominal_state[i] += K_yaw[i] * innovation_yaw;
        }
        nominal_state[4] = normalizeAngle(nominal_state[4]); // Chuẩn hóa lại góc sau khi cập nhật

        // Cập nhật hiệp phương sai: P = (I - K * H) * P
        double I_KH[STATE_DIM * STATE_DIM];
        for (int i = 0; i < STATE_DIM; i++)
        {
          for (int j = 0; j < STATE_DIM; j++)
          {
            // H_yaw chỉ có phần tử thứ 4 là khác 0
            double kh_ij = K_yaw[i] * ((j == 4) ? 1.0 : 0.0);
            I_KH[i * STATE_DIM + j] = ((i == j) ? 1.0 : 0.0) - kh_ij;
          }
        }

        double P_new[STATE_DIM * STATE_DIM] = {0};
        for (int i = 0; i < STATE_DIM; i++)
        {
          for (int j = 0; j < STATE_DIM; j++)
          {
            for (int k = 0; k < STATE_DIM; k++)
            {
              P_new[i * STATE_DIM + j] += I_KH[i * STATE_DIM + k] * P[k * STATE_DIM + j];
            }
          }
        }
        memcpy(P, P_new, sizeof(P_new));
      }
      // Sau mỗi bước, in ra state
      // Serial.print("x=");
      // Serial.print(nominal_state[0], 3);
      // Serial.print("  y=");
      // Serial.print(nominal_state[1], 3);
      // Serial.print("  yaw=");
      // Serial.println(nominal_state[4] * 180.0 / M_PI, 2);
      xSemaphoreGive(stateMutex);
    }
    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
}

void updateTask(void *pvParameters)
{
  const TickType_t xFrequency = 10 / portTICK_PERIOD_MS; // 100Hz
  TickType_t xLastWakeTime = xTaskGetTickCount();
  uint8_t rtksta = 0;
  uint32_t hAcc = 0;
  PointECEF tempECEF;
  PointECEF tempECEF_2;
  while (1)
  {
    bool NAVHPPOSECEF = false;
    bool PVT = false;

    double z[2] = {0, 0}; // Vẫn là [current_east, current_north]
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
      // double current_east_2, current_north_2, current_up_2;
      // ecefToEnu(tempECEF_2, current_east_2, current_north_2, current_up_2);
      z[0] = current_east;
      z[1] = current_north;

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
      hAcc = myGNSS.getHorizontalAccEst(); // mm
      NAVHPPOSECEF = true;                 // Đánh dấu đã có dữ liệu hợp lệ
      // PVT = true;
    }

    // if (myGNSS.getPVT() && (myGNSS.getInvalidLlh() == false))
    // {
    //   rtksta = myGNSS.getCarrierSolutionType();
    //   hAcc = myGNSS.getHorizontalAccEst(); // mm
    //   PVT = true;
    // }

    if (NAVHPPOSECEF)
    {
      if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
      {
        Serial.print(z[0]);
        Serial.print(",");
        Serial.print(z[1]);
        Serial.print(",  ");

        // Serial.print("Pred");
        Serial.print(nominal_state[0], 4);
        Serial.print(",");
        Serial.print(nominal_state[1], 4);
        Serial.print(" , ");
        Serial.print(degrees(nominal_state[4]), 4);
        Serial.print(" , ");

        // =================== BƯỚC 3: TÍNH TOÁN HIỆU CHỈNH LỖI ===================

        // Ma trận H cho trạng thái lỗi vẫn như cũ
        static double H[2 * STATE_DIM] = {0};
        H[0 * STATE_DIM + 0] = 1.0; // d(z_e)/d(δp_e)
        H[1 * STATE_DIM + 1] = 1.0; // d(z_n)/d(δp_n)

        // =================== BẮT ĐẦU ADAPTIVE R ===================
        double sigma = (double)hAcc / 1000.0; // hAcc từ mm sang m
        static double R_local[4] = {0};
        Serial.print(hAcc);
        Serial.print(", ");
        R_local[0] = sigma * sigma; // phương sai East
        R_local[3] = sigma * sigma; // phương sai North
        // Serial.print(R_local[0], 5);
        // Serial.print(", ");
        // =================== KẾT THÚC ADAPTIVE R ===================

        // =================== TÍNH TOÁN TỐI ƯU ===================
        // --- BƯỚC 1: TÍNH S = H*P*H^T + R ---
        static double S[4];
        S[0] = P[0 * STATE_DIM + 0] + R_local[0]; // S(0,0)
        S[1] = P[0 * STATE_DIM + 1] + R_local[1]; // S(0,1)
        S[2] = P[1 * STATE_DIM + 0] + R_local[2]; // S(1,0)
        S[3] = P[1 * STATE_DIM + 1] + R_local[3]; // S(1,1)

        // --- BƯỚC 2: TÍNH S_inv và K = P*H^T*S_inv ---
        double detS = S[0] * S[3] - S[1] * S[2];
        if (fabs(detS) < 1e-9)
        {
          //*********DEBUG***********//
          // Serial.print("detS:");
          // Serial.println(detS, 10);
          //*********DEBUG***********//

          Serial.println("Warning: det(S) is too small. Skipping EKF update.");
          xSemaphoreGive(stateMutex);
        }
        else
        {
          // 3a. Tính S_inv (kích thước 2x2)
          static double S_inv[4];
          S_inv[0] = S[3] / detS;
          S_inv[1] = -S[1] / detS;
          S_inv[2] = -S[2] / detS;
          S_inv[3] = S[0] / detS;

          static double K[STATE_DIM * 2];
          for (int i = 0; i < STATE_DIM; i++)
          {
            double p_i0 = P[i * STATE_DIM + 0];
            double p_i1 = P[i * STATE_DIM + 1];
            K[i * 2 + 0] = p_i0 * S_inv[0] + p_i1 * S_inv[2];
            K[i * 2 + 1] = p_i0 * S_inv[1] + p_i1 * S_inv[3];
          }

          // --- BƯỚC 3: CẬP NHẬT TRẠNG THÁI ---
          double innovation[2] = {z[0] - nominal_state[0], z[1] - nominal_state[1]};

          // // d_squared = innovation' * S_inv * innovation
          // double d_squared = S_inv[0] * innovation[0] * innovation[0] +
          //                    (S_inv[1] + S_inv[2]) * innovation[0] * innovation[1] +
          //                    S_inv[3] * innovation[1] * innovation[1];

          // const double CHI_SQUARED_THRESHOLD = 5.991; // Ngưỡng 95% cho 2 bậc tự do

          // if (d_squared > CHI_SQUARED_THRESHOLD)
          // {
          //   Serial.println("Outlier detected by Chi-Squared test! Skipping update.");
          //   xSemaphoreGive(stateMutex);
          //   // vTaskDelayUntil(&xLastWakeTime, xFrequency);
          //   continue;
          // }
          // Tính toán véc-tơ hiệu chỉnh lỗi (δx = K * innovation)
          double error_state_correction[STATE_DIM] = {0};
          // error_state_correction = K * innovation
          for (int i = 0; i < STATE_DIM; i++)
          {
            error_state_correction[i] = K[i * 2 + 0] * innovation[0] + K[i * 2 + 1] * innovation[1];
          }

          // =================== BƯỚC 4: TIÊM LỖI VÀ RESET ===================

          // "Tiêm" (Inject) sự hiệu chỉnh lỗi vào trạng thái danh nghĩa
          for (int i = 0; i < STATE_DIM; i++)
          {
            nominal_state[i] += error_state_correction[i];
          }
          // Chuẩn hóa lại góc sau khi hiệu chỉnh
          nominal_state[4] = normalizeAngle(nominal_state[4]);

          // Đây là cách tính (I - K*H)*P hiệu quả nhất
          static double P_new[STATE_DIM * STATE_DIM];
          for (int i = 0; i < STATE_DIM; i++)
          {
            for (int j = 0; j < STATE_DIM; j++)
            {
              double kh_ij = K[i * 2 + 0] * P[0 * STATE_DIM + j] + K[i * 2 + 1] * P[1 * STATE_DIM + j];
              P_new[i * STATE_DIM + j] = P[i * STATE_DIM + j] - kh_ij;
            }
          }
          memcpy(P, P_new, sizeof(P_new));

          // Trạng thái lỗi đã được "tiêu thụ" và được reset về 0 cho vòng lặp tiếp theo.
          // Hiệp phương sai P đã được cập nhật để phản ánh sự không chắc chắn mới.
          Serial.print(rtksta);
          Serial.print(", ");
          Serial.print(nominal_state[0], 4);
          Serial.print(",");
          Serial.print(nominal_state[1], 4);
          Serial.print(",");
          Serial.print(nominal_state[4], 4);
          Serial.print(",");
          Serial.println(nominal_state[4] * 180 / M_PI, 4);
          xSemaphoreGive(stateMutex);
        }
      }
    }
    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
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

void setup()
{
  setCpuFrequencyMhz(240);
  delay(500);

  Serial.begin(115200);
  delay(500);

  Serial1.begin(115200, SERIAL_8N1, RXD2, TXD2);
  delay(500);

  WitInit(WIT_PROTOCOL_MODBUS, 0x50);
  WitSerialWriteRegister(SensorUartSend);
  WitRegisterCallBack(CopeSensorData);
  WitDelayMsRegister(Delayms);
  // WitSetHeadingtoZero();

  // WiFi.mode(WIFI_STA);
  // delay(500);

  // Init ESP-NOW
  // if (esp_now_init() != ESP_OK)
  // {
  //   Serial.println("Error initializing ESP-NOW");
  //   return;
  // }
  // Serial.println("222222");
  // memcpy(peerInfo.peer_addr, broadcastAddress, 6);
  // peerInfo.channel = 0;
  // peerInfo.encrypt = false;
  // if (esp_now_add_peer(&peerInfo) != ESP_OK)
  // {
  //   Serial.println("Failed to add peer");
  //   return;
  // }

  Wire.begin(21, 22);
  delay(500);

  if (myGNSS.begin(Wire, 0x42) == false)
  {
    Serial.println(F("u-blox GNSS not detected. Freezing."));
    while (1)
      ;
  }

  stateMutex = xSemaphoreCreateMutex();

  myGNSS.setI2COutput(COM_TYPE_UBX);

  // sensor.setCalibrationCommand(SET_HEADING_ANGLE_TO_ZERO);
  // Chờ RTK đạt mức Fixed (fix type = 3)
  Serial.println("Waiting for RTK Fixed solution...");
  while (true)
  {
    // Serial.println(sensor.getAttitudeValues().yaw);
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

      GetRTK(initialLLH);
      // Thiết lập gốc ENU
      CalculateEnuMatrix(initialLLH);
    }
    if (myGNSS.getNAVHPPOSECEF())
    {
      enu_origin_ecef.x = ((double)myGNSS.getHighResECEFX() / 100.0) + ((double)myGNSS.getHighResECEFXHp() / 10000.0);
      enu_origin_ecef.y = ((double)myGNSS.getHighResECEFY() / 100.0) + ((double)myGNSS.getHighResECEFYHp() / 10000.0);
      enu_origin_ecef.z = ((double)myGNSS.getHighResECEFZ() / 100.0) + ((double)myGNSS.getHighResECEFZHp() / 10000.0);
      Serial.print("ENU Origin ECEF: ");
      Serial.print(enu_origin_ecef.x, 4);
      Serial.print(", ");
      Serial.print(enu_origin_ecef.y, 4);
      Serial.print(", ");
      Serial.println(enu_origin_ecef.z, 4);
    }
    if (fixType >= 3 && carrierSolution >= 1)
    { // RTK Fixed
      Serial.println("RTK Fixed acquired.");
      break;
    }

    delay(100);
  }

  nominal_state[0] = 0.0; // x ban đầu
  nominal_state[1] = 0.0; // y ban đầu
  // prev_rtk_pos[0] = initialECEF.x;
  // prev_rtk_pos[1] = initialECEF.y;
  Serial.print("Initial RTK: x=");
  Serial.print(nominal_state[0]);
  Serial.print(", y=");
  Serial.print(nominal_state[1]);

  nominal_state[2] = 0.0; // v_x = 0
  nominal_state[3] = 0.0; // v_y = 0
  // Lấy heading ban đầu từ WTGAHRS3
  // AttitudeData initialHeading = sensor.getAttitudeValues(); // Giả sử hàm này trả về độ (-180 đến 180)
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
    if (s_cDataUpdate & ANGLE_UPDATE)
    {
      nominal_state[4] = fAngle[2];
      s_cDataUpdate &= ~ANGLE_UPDATE;
    }
    s_cDataUpdate = 0;
  }
  else
  {
    nominal_state[4] = 0; // Mặc định 0 nếu không đọc được
    Serial.println("Warning: Initial heading not available, set to 0.");
  }

  // if (initialHeading.isDataValid != 0)
  // {                                                       // Kiểm tra dữ liệu hợp lệ
  //   nominal_state[4] = initialHeading.yaw * M_PI / 180.0; // Chuyển sang radian
  // }

  Serial.print("Initial Heading (rad): ");
  Serial.println(nominal_state[4]);
  nominal_state[5] = 0; // b_ax = 0
  nominal_state[6] = 0; // b_ay = 0
  nominal_state[7] = 0; // b_gz = 0

  if (stateMutex == NULL)
  {
    Serial.println("Error: Failed to create mutex!");
    while (1)
      ; // Treo hệ thống
  }
  double yaw_measurement_std_dev_rad = 0.1 * M_PI / 180.0;
  R_yaw = yaw_measurement_std_dev_rad * yaw_measurement_std_dev_rad;

  // =================== BỘ THAM SỐ ĐỀ XUẤT ===================
  // Khởi tạo P - Mức độ không chắc chắn BAN ĐẦU
  P[0 * STATE_DIM + 0] = 0.1 * 0.1;                                   // Vị trí East (độ lệch chuẩn 10cm)
  P[1 * STATE_DIM + 1] = 0.1 * 0.1;                                   // Vị trí North (độ lệch chuẩn 10cm)
  P[2 * STATE_DIM + 2] = 0.1 * 0.1;                                   // Vận tốc East (độ lệch chuẩn 10cm/s)
  P[3 * STATE_DIM + 3] = 0.1 * 0.1;                                   // Vận tốc North (độ lệch chuẩn 10cm/s)
  P[4 * STATE_DIM + 4] = (5.0 * M_PI / 180.0) * (5.0 * M_PI / 180.0); // Heading (độ lệch chuẩn 5 độ)
  P[5 * STATE_DIM + 5] = 1.0;                                         // Rất không chắc chắn về bias gia tốc ban đầu X
  P[6 * STATE_DIM + 6] = 1.0;                                         // Rất không chắc chắn về bias gia tốc ban đầu Y
  P[7 * STATE_DIM + 7] = (0.5 * M_PI / 180.0) * (0.5 * M_PI / 180.0); // Không chắc chắn về bias gyro ban đầu (0.5 deg/s)

  // Khởi tạo Q - Mức độ nhiễu của quá trình THEO THỜI GIAN
  double vel_noise_std_dev = 0.05;             // Giả định vận tốc có thể nhiễu ngẫu nhiên 5 cm/s
  double accel_bias_noise_std_dev = 0.0026999; // Giả định bias gia tốc RẤT ỔN ĐỊNH, ít thay đổi
  double gyro_bias_noise_std_dev = 0.01994;   // Giả định bias gyro CỰC KỲ ỔN ĐỊNH

  Q[2 * STATE_DIM + 2] = vel_noise_std_dev * vel_noise_std_dev;
  Q[3 * STATE_DIM + 3] = vel_noise_std_dev * vel_noise_std_dev;
  Q[4 * STATE_DIM + 4] = (0.01 * M_PI / 180.0) * (0.01 * M_PI / 180.0); // Nhiễu đo của Gyro (từ datasheet)
  Q[5 * STATE_DIM + 5] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;
  Q[6 * STATE_DIM + 6] = accel_bias_noise_std_dev * accel_bias_noise_std_dev;
  Q[7 * STATE_DIM + 7] = gyro_bias_noise_std_dev * gyro_bias_noise_std_dev;

  // Tạo tác vụ cho hai core

  xTaskCreatePinnedToCore(predictionTask, "Prediction", 4096, NULL, 1, NULL, 0); // Core 0
  xTaskCreatePinnedToCore(updateTask, "Update", 4096, NULL, 1, NULL, 1);         // Core 1

  Serial.println("System initialized with RTK Fixed position and WTGAHRS3 heading.");
  delay(2000);
}

void loop() {}