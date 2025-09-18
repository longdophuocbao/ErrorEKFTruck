// ES-EKF
#include <Wire.h>
#include <SparkFun_u-blox_GNSS_v3.h>
#include <ModbusMaster.h>
#include <WTGAHRS3_485.h>
#include <esp_now.h>
#include <WiFi.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>

#define RXD2 16
#define TXD2 17
#define SENSOR_SLAVE_ID 0x50

// --- CẤU HÌNH BỘ LỌC ---
const int STATE_DIM = 8; // Kích thước của cả Nominal State và Error State

// --- BIẾN TOÀN CỤC ---
SFE_UBLOX_GNSS myGNSS;

esp_now_peer_info_t peerInfo;
uint8_t broadcastAddress[] = {0x24, 0xD7, 0xEB, 0x13, 0x74, 0xAC};

ModbusMaster node;
HardwareSerial RS485(1);
WTGAHRS3_485 sensor(node);
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

// Trạng thái danh nghĩa (Nominal State) - Ước tính chính của hệ thống
double nominal_state[STATE_DIM] = {0}; // [p_e, p_n, v_e, v_n, theta, b_ax, b_ay, b_gz]

// Ma trận hiệp phương sai của TRẠNG THÁI LỖI (Error State Covariance)
double P[STATE_DIM * STATE_DIM] = {0};

// Ma trận nhiễu quá trình cho TRẠNG THÁI LỖI
double Q[STATE_DIM * STATE_DIM] = {0};

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
  long latitude = myGNSS.getHighResLatitude();
  int8_t latitudeHp = myGNSS.getHighResLatitudeHp();
  _CurrentPos.latitude = ((double)latitude / 10000000.0) + ((double)latitudeHp / 1000000000.0);

  long longitude = myGNSS.getHighResLongitude();
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

  east = ecef_to_enu_matrix[0][0] * dx + ecef_to_enu_matrix[0][1] * dy + ecef_to_enu_matrix[0][2] * dz;
  north = ecef_to_enu_matrix[1][0] * dx + ecef_to_enu_matrix[1][1] * dy + ecef_to_enu_matrix[1][2] * dz;
  up = ecef_to_enu_matrix[2][0] * dx + ecef_to_enu_matrix[2][1] * dy + ecef_to_enu_matrix[2][2] * dz;
}

/**
 * @brief Thiết lập gốc ENU và ma trận chuyển đổi từ ECEF sang ENU.
 * @param origin_llh Tọa độ gốc ENU ở dạng LLH.
 * @param origin_ecef Tọa độ gốc ENU ở dạng ECEF.
 */
void setupEnuOrigin(const PointLLH &origin_llh, const PointECEF &origin_ecef)
{
  enu_origin_ecef = origin_ecef;

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
}

/**
 * @brief Chuẩn hóa góc về khoảng [-π, π].
 * @param angle Góc đầu vào (radian).
 * @return Góc đã chuẩn hóa (radian).
 */
double normalizeAngle(double angle)
{
  // Dùng fmod để đưa góc về khoảng [-2π, 2π]
  angle = fmod(angle + M_PI, 2.0 * M_PI);
  if (angle < 0)
  {
    angle += 2.0 * M_PI;
  }
  return angle - M_PI;
}
// --- KẾT THÚC KHAI BÁO ---

void predictionTask(void *pvParameters)
{
  const TickType_t xFrequency = 10 / portTICK_PERIOD_MS; // 100Hz
  TickType_t xLastWakeTime = xTaskGetTickCount();
  while (1)
  {
    // Đọc IMU và chuyển đổi gyro_z sang rad/s
    AccelerationData Accel = sensor.getAccelerationData();
    AngularVelocityData Gyro = sensor.getAngularVelocityData();
    if (Accel.isDataValid && Gyro.isDataValid)
    {
      accel[0] = Accel.accelX;
      accel[1] = Accel.accelY;
      gyro_z = Gyro.angularVelZ * M_PI / 180.0;
    }
    else
    {
      continue;
    }
    double dt = 0.01;

    if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
    {
      // =================== BƯỚC 1: CẬP NHẬT TRẠNG THÁI DANH NGHĨA (NOMINAL STATE) ===================
      // Sử dụng RK4 để tích hợp IMU vào nominal_state, làm cho nó bị trôi đi một cách tự nhiên.
      // (Code RK4 giữ nguyên như phiên bản trước, nhưng tác động lên nominal_state)
      {
        double ax_raw = accel[0] - nominal_state[5];
        double ay_raw = accel[1] - nominal_state[6];
        double gyro_z_raw = gyro_z - nominal_state[7];

        double ax_std_body = -ay_raw;
        double ay_std_body = ax_raw;

        // 3. Lấy góc heading (yaw) từ state
        double theta = nominal_state[4];
        double cos_th = cos(theta);
        double sin_th = sin(theta);

        // 4. Xoay gia tốc từ Body tiêu chuẩn sang ENU(East, North)
        // Chú ý công thức đúng: East ~ Y, North ~ X trong mặt phẳng toán học
        // double ax_nav_east = ay_std_body * sin_th + ax_std_body * cos_th;
        // double ay_nav_north = ay_std_body * cos_th - ax_std_body * sin_th;
        double ax_nav_east = ax_std_body * cos_th - ay_std_body * sin_th;
        double ay_nav_north = ax_std_body * sin_th + ay_std_body * cos_th;

        double k1[STATE_DIM], k2[STATE_DIM], k3[STATE_DIM], k4[STATE_DIM];
        double temp_state[STATE_DIM];

        // -- Bước 1: Tính k1 dựa trên trạng thái hiện tại (state)
        k1[0] = nominal_state[2];
        k1[1] = nominal_state[3];
        k1[2] = ax_nav_east;
        k1[3] = ay_nav_north;
        k1[4] = gyro_z_raw;
        k1[5] = 0; // Đạo hàm của bias là 0
        k1[6] = 0;
        k1[7] = 0;

        // -- Bước 2: Tính k2 dựa trên trạng thái ở điểm giữa (sử dụng k1)
        for (int i = 0; i < STATE_DIM; ++i)
          temp_state[i] = nominal_state[i] + 0.5 * dt * k1[i];
        cos_th = cos(temp_state[4]);
        sin_th = sin(temp_state[4]);
        ax_nav_east = ax_std_body * cos_th - ay_std_body * sin_th;
        ay_nav_north = ax_std_body * sin_th + ay_std_body * cos_th;
        k2[0] = temp_state[2];
        k2[1] = temp_state[3];
        k2[2] = ax_nav_east;
        k2[3] = ay_nav_north;
        k2[4] = gyro_z_raw;
        for (int i = 5; i < STATE_DIM; ++i)
          k2[i] = 0;

        // -- Bước 3: Tính k3 dựa trên trạng thái ở điểm giữa (sử dụng k2)
        for (int i = 0; i < STATE_DIM; ++i)
          temp_state[i] = nominal_state[i] + 0.5 * dt * k2[i];
        cos_th = cos(temp_state[4]);
        sin_th = sin(temp_state[4]);
        ax_nav_east = ax_std_body * cos_th - ay_std_body * sin_th;
        ay_nav_north = ax_std_body * sin_th + ay_std_body * cos_th;
        k3[0] = temp_state[2];
        k3[1] = temp_state[3];
        k3[2] = ax_nav_east;
        k3[3] = ay_nav_north;
        k3[4] = gyro_z_raw;
        for (int i = 5; i < STATE_DIM; ++i)
          k3[i] = 0;

        // -- Bước 4: Tính k4 dựa trên trạng thái ở điểm cuối (sử dụng k3)
        for (int i = 0; i < STATE_DIM; ++i)
          temp_state[i] = nominal_state[i] + dt * k3[i];
        cos_th = cos(temp_state[4]);
        sin_th = sin(temp_state[4]);
        ax_nav_east = ax_std_body * cos_th - ay_std_body * sin_th;
        ay_nav_north = ax_std_body * sin_th + ay_std_body * cos_th;
        k4[0] = temp_state[2];
        k4[1] = temp_state[3];
        k4[2] = ax_nav_east;
        k4[3] = ay_nav_north;
        k4[4] = gyro_z_raw;
        for (int i = 5; i < STATE_DIM; ++i)
          k4[i] = 0;

        // Bước 4: Cập nhật trạng thái
        for (int i = 0; i < STATE_DIM; i++)
        {
          nominal_state[i] += dt * k2[i];
        }

        // Chuẩn hóa lại góc theta
        nominal_state[4] = normalizeAngle(nominal_state[4]);

        xSemaphoreGive(stateMutex);
      }

      // =================== BƯỚC 2: CẬP NHẬT HIỆP PHƯƠNG SAI LỖI (ERROR STATE COVARIANCE) ===================
      // P_k = F * P_{k-1} * F^T + Q
      // Chúng ta cần tính ma trận chuyển trạng thái lỗi F

      double F[STATE_DIM * STATE_DIM] = {0};
      double theta = nominal_state[4];
      double cos_th = cos(theta);
      double sin_th = sin(theta);

      // Gia tốc đã hiệu chỉnh bias trong hệ body
      double corrected_ax_raw = accel[0] - nominal_state[5];
      double corrected_ay_raw = accel[1] - nominal_state[6];
      double ax_std_body = -corrected_ay_raw;
      double ay_std_body = corrected_ax_raw;

      // Xây dựng ma trận F ≈ I + A*dt
      for (int i = 0; i < STATE_DIM; ++i)
        F[i * STATE_DIM + i] = 1.0; // Bắt đầu với ma trận đơn vị I

      F[0 * STATE_DIM + 2] = dt; // d(p_e)/d(v_e)
      F[1 * STATE_DIM + 3] = dt; // d(p_n)/d(v_n)

      // d(v_e)/d(theta)
      F[2 * STATE_DIM + 4] = (-ax_std_body * sin_th - ay_std_body * cos_th) * dt;
      // d(v_e)/d(b_ax)
      F[2 * STATE_DIM + 5] = (sin_th)*dt; // d(v_e)/d(b_ay) = d(v_e)/d(ax_std) * d(ax_std)/d(bay) = -cos(th) * -1
      // d(v_e)/d(b_ay)
      F[2 * STATE_DIM + 6] = (-cos_th) * dt;

      // d(v_n)/d(theta)
      F[3 * STATE_DIM + 4] = (ax_std_body * cos_th - ay_std_body * sin_th) * dt;
      // d(v_n)/d(b_ax)
      F[3 * STATE_DIM + 5] = (-cos_th) * dt;
      // d(v_n)/d(b_ay)
      F[3 * STATE_DIM + 6] = (-sin_th) * dt;

      F[4 * STATE_DIM + 7] = -dt; // d(theta)/d(b_gz)

      // Cập nhật P: P = F * P * F^T
      double FP[STATE_DIM * STATE_DIM] = {0};
      // FP = F * P
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
      // P = FP * F^T + Q
      for (int i = 0; i < STATE_DIM; i++)
      {
        for (int j = 0; j < STATE_DIM; j++)
        {
          double temp = 0;
          for (int k = 0; k < STATE_DIM; k++)
          {
            temp += FP[i * STATE_DIM + k] * F[j * STATE_DIM + k]; // F^T(k,j) = F(j,k)
          }
          P[i * STATE_DIM + j] = temp + Q[i * STATE_DIM + j];
        }
      }

      xSemaphoreGive(stateMutex);
    }
    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
}

void updateTask(void *pvParameters)
{
  const TickType_t xFrequency = 1000 / portTICK_PERIOD_MS; // 1Hz
  TickType_t xLastWakeTime = xTaskGetTickCount();
  uint8_t rtksta = 0;
  while (1)
  {
    bool gnss_data_valid = false;
    double z[2] = {0, 0}; // Vẫn là [current_east, current_north]

    if (myGNSS.getPVT())
    {
      rtksta = myGNSS.getCarrierSolutionType();
      PointLLH tempLLH;
      PointECEF tempECEF;

      // Lấy tọa độ và chuyển sang hệ ECEF
      GetRTK(tempLLH);
      llhToEcef(tempLLH, tempECEF);

      // BƯỚC 2: Chuyển đổi tọa độ ECEF vừa nhận được sang ENU
      double current_east, current_north, current_up;
      ecefToEnu(tempECEF, current_east, current_north, current_up);

      z[0] = current_east;
      z[1] = current_north;

      //*********DEBUG***********//
      Serial.print(current_east);
      Serial.print(", ");
      Serial.print(current_north);
      Serial.print(", ");
      //*********DEBUG***********//

      gnss_data_valid = true; // Đánh dấu đã có dữ liệu hợp lệ
    }

    if (gnss_data_valid)
    {
      if (xSemaphoreTake(stateMutex, (TickType_t)10) == pdTRUE)
      {
        Serial.print("Prediction: state[0]=");
        Serial.print(nominal_state[0], 6);
        Serial.print(", state[1]=");
        Serial.print(nominal_state[1], 6);
        Serial.print(" , ");

        // =================== BƯỚC 3: TÍNH TOÁN HIỆU CHỈNH LỖI ===================

        // Innovation: Sai số giữa phép đo GPS và trạng thái danh nghĩa
        double innovation[2] = {z[0] - nominal_state[0], z[1] - nominal_state[1]};

        // Ma trận H cho trạng thái lỗi vẫn như cũ
        double H[2 * STATE_DIM] = {0};
        H[0 * STATE_DIM + 0] = 1.0; // d(z_e)/d(δp_e)
        H[1 * STATE_DIM + 1] = 1.0; // d(z_n)/d(δp_n)

        // =================== BẮT ĐẦU ADAPTIVE R ===================
        double R_local[2 * 2]; // Khai báo R_local là biến cục bộ trong task

        switch (rtksta)
        {
        case 2: // RTK Fixed - Độ chính xác cao nhất
          // Giả định độ lệch chuẩn ~2cm -> phương sai = 0.02*0.02 = 0.0004
          R_local[0] = 0.0004;
          R_local[1] = 0;
          R_local[2] = 0;
          R_local[3] = 0.0004;
          break;

        case 1: // RTK Float - Độ chính xác trung bình
          // Giả định độ lệch chuẩn ~30cm -> phương sai = 0.3*0.3 = 0.09
          R_local[0] = 0.09;
          R_local[1] = 0;
          R_local[2] = 0;
          R_local[3] = 0.09;
          break;

        default: // Trường hợp khác (No RTK, Standard GPS) - Độ chính xác thấp
          // Giả định độ lệch chuẩn ~1.5m -> phương sai = 1.5*1.5 = 2.25
          R_local[0] = 2.25;
          R_local[1] = 0;
          R_local[2] = 0;
          R_local[3] = 2.25;
          break;
        }
        // =================== KẾT THÚC ADAPTIVE R ===================

        // =================== TÍNH TOÁN TỐI ƯU ===================
        // --- BƯỚC 1: TÍNH S = H*P*H^T + R ---
        double S[4];
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
          double S_inv[4];
          S_inv[0] = S[3] / detS;
          S_inv[1] = -S[1] / detS;
          S_inv[2] = -S[2] / detS;
          S_inv[3] = S[0] / detS;

          double K[STATE_DIM * 2];
          for (int i = 0; i < STATE_DIM; i++)
          {
            double p_i0 = P[i * STATE_DIM + 0];
            double p_i1 = P[i * STATE_DIM + 1];
            K[i * 2 + 0] = p_i0 * S_inv[0] + p_i1 * S_inv[2];
            K[i * 2 + 1] = p_i0 * S_inv[1] + p_i1 * S_inv[3];
          }

          // --- BƯỚC 3: CẬP NHẬT TRẠNG THÁI ---
          double innovation[2] = {z[0] - nominal_state[0], z[1] - nominal_state[1]};
          for (int i = 0; i < STATE_DIM; i++)
          {
            nominal_state[i] += K[i * 2 + 0] * innovation[0] + K[i * 2 + 1] * innovation[1];
          }

          // Tính toán véc-tơ hiệu chỉnh lỗi tối ưu (δx)
          double error_state_correction[STATE_DIM] = {0};
          // error_state_correction = K * innovation
          for (int i = 0; i < STATE_DIM; i++)
          {
            error_state_correction[i] = K[i * 2 + 0] * innovation[0] + K[i * 2 + 1] * innovation[1];
          }

          // --- BƯỚC 4: CẬP NHẬT HIỆP PHƯƠNG SAI P = (I - K*H)*P ---
          // Đây là cách tính (I - K*H)*P hiệu quả nhất
          double P_new[STATE_DIM * STATE_DIM];
          for (int i = 0; i < STATE_DIM; i++)
          {
            for (int j = 0; j < STATE_DIM; j++)
            {
              double kh_ij = K[i * 2 + 0] * P[0 * STATE_DIM + j] + K[i * 2 + 1] * P[1 * STATE_DIM + j];
              P_new[i * STATE_DIM + j] = P[i * STATE_DIM + j] - kh_ij;
            }
          }
          memcpy(P, P_new, sizeof(P_new));

          // =================== BƯỚC 4: TIÊM LỖI VÀ RESET ===================

          // "Tiêm" (Inject) sự hiệu chỉnh lỗi vào trạng thái danh nghĩa
          for (int i = 0; i < STATE_DIM; i++)
          {
            nominal_state[i] += error_state_correction[i];
          }
          // Chuẩn hóa lại góc sau khi hiệu chỉnh
          nominal_state[4] = normalizeAngle(nominal_state[4]);

          // Trạng thái lỗi đã được "tiêu thụ" và được reset về 0 cho vòng lặp tiếp theo.
          // Hiệp phương sai P đã được cập nhật để phản ánh sự không chắc chắn mới.
          Serial.print(rtksta);
          Serial.print(",");
          Serial.print(nominal_state[0], 4);
          Serial.print(",");
          Serial.print(nominal_state[1], 4);
          Serial.print(",");
          Serial.println(nominal_state[4], 4);
          xSemaphoreGive(stateMutex);
        }
      }
    }
    vTaskDelayUntil(&xLastWakeTime, xFrequency);
  }
}

void setup()
{
  setCpuFrequencyMhz(240);
  Serial.begin(115200);
  delay(500);
  Serial.println("SparkFun u-blox Example");
  RS485.begin(115200, SERIAL_8N1, RXD2, TXD2);
  delay(500);
  node.begin(SENSOR_SLAVE_ID, RS485);
  delay(500);
  WiFi.mode(WIFI_STA);
  delay(500);

  // Init ESP-NOW
  if (esp_now_init() != ESP_OK)
  {
    Serial.println("Error initializing ESP-NOW");
    return;
  }
  memcpy(peerInfo.peer_addr, broadcastAddress, 6);
  peerInfo.channel = 0;
  peerInfo.encrypt = false;
  if (esp_now_add_peer(&peerInfo) != ESP_OK)
  {
    Serial.println("Failed to add peer");
    return;
  }

  Wire.begin(21, 22);
  delay(500);
  if (myGNSS.begin(Wire, 0x42) == false)
  {
    Serial.println(F("u-blox GNSS not detected. Freezing."));
    while (1)
      ;
  }

  stateMutex = xSemaphoreCreateMutex();

  // Chờ RTK đạt mức Fixed (fix type = 3)
  Serial.println("Waiting for RTK Fixed solution...");
  while (true)
  {
    if (myGNSS.getPVT())
    {
      uint8_t fixType = myGNSS.getFixType();
      uint8_t carrierSolution = myGNSS.getCarrierSolutionType();
      Serial.print("Fix Type: ");
      Serial.print(fixType);
      Serial.print(", Carrier Solution: ");
      Serial.println(carrierSolution);
      if (fixType >= 3 && carrierSolution >= 1)
      { // RTK Fixed
        Serial.println("RTK Fixed acquired.");
        break;
      }
    }
    delay(100);
  }

  // Lấy vị trí ban đầu từ RTK
  PointLLH initialLLH;
  PointECEF initialECEF;
  GetRTK(initialLLH);
  llhToEcef(initialLLH, initialECEF);
  // Thiết lập gốc ENU
  setupEnuOrigin(initialLLH, initialECEF);

  // state[0] = initialECEF.x; // x ban đầu
  // state[1] = initialECEF.y; // y ban đầu
  // prev_rtk_pos[0] = initialECEF.x;
  // prev_rtk_pos[1] = initialECEF.y;
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
  AttitudeData initialHeading = sensor.getAttitudeValues(); // Giả sử hàm này trả về độ (-180 đến 180)
  if (initialHeading.isDataValid != 0)
  {                                                       // Kiểm tra dữ liệu hợp lệ
    nominal_state[4] = initialHeading.yaw * M_PI / 180.0; // Chuyển sang radian
  }
  else
  {
    nominal_state[4] = 0; // Mặc định 0 nếu không đọc được
    Serial.println("Warning: Initial heading not available, set to 0.");
  }
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
  double vel_noise_std_dev = 0.05;        // Giả định vận tốc có thể nhiễu ngẫu nhiên 5 cm/s
  double accel_bias_noise_std_dev = 0.01; // Giả định bias gia tốc RẤT ỔN ĐỊNH, ít thay đổi
  double gyro_bias_noise_std_dev = 0.001; // Giả định bias gyro CỰC KỲ ỔN ĐỊNH

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