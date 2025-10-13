import os
import subprocess
from SCons.Script import DefaultEnvironment

# Lấy môi trường build của PlatformIO
env = DefaultEnvironment()

# Lấy đường dẫn tới thư mục chứa các thư viện của dự án
lib_deps_dir = env.subst("$PROJECT_LIBDEPS_DIR")
pio_env = env.subst("$PIOENV")
nanopb_lib_name = "Nanopb"

# Tạo đường dẫn đầy đủ tới thư mục thư viện Nanopb
nanopb_lib_dir = os.path.join(lib_deps_dir, pio_env, nanopb_lib_name)

# Kiểm tra xem thư mục thư viện có thực sự tồn tại không
if not os.path.isdir(nanopb_lib_dir):
    print(f"Error: Could not find the Nanopb library directory.")
    print(f"Looked for: {nanopb_lib_dir}")
    env.Exit(1)

# Đường dẫn tới script nanopb_generator.py
nanopb_generator = os.path.join(nanopb_lib_dir, "generator", "nanopb_generator.py")
# --- PHẦN SỬA LỖI ---
# Đường dẫn tới thư mục chứa file nanopb.proto (để import)
nanopb_proto_include_dir = os.path.join(nanopb_lib_dir, "generator", "proto")
# --- KẾT THÚC PHẦN SỬA LỖI ---


# Đường dẫn tới thư mục chứa các file .proto của bạn
PROTO_DIR = os.path.join(env.subst("$PROJECT_DIR"), "proto")

# Đường dẫn tới thư mục sẽ chứa các file .c và .h được tạo ra
GENERATED_DIR = os.path.join(env.subst("$PROJECT_SRC_DIR"), "generated")

# Tạo thư mục generated nếu nó chưa tồn tại
if not os.path.exists(GENERATED_DIR):
    os.makedirs(GENERATED_DIR)

# Lặp qua tất cả các file .proto trong thư mục proto
for filename in os.listdir(PROTO_DIR):
    if filename.endswith(".proto"):
        proto_file_path = os.path.join(PROTO_DIR, filename)
        print(f"Generating C/C++ files for {filename}...")
        
        # --- PHẦN SỬA LỖI ---
        # Cập nhật lại lệnh để thêm các đường dẫn -I
        command = [
            env.subst("$PYTHONEXE"),
            nanopb_generator,
            # Chỉ định thư mục đầu ra cho file .pb.c và .pb.h
            f"-D",
            GENERATED_DIR,
            # Thêm thư mục chứa file nanopb.proto vào đường dẫn tìm kiếm
            f"-I{nanopb_proto_include_dir}",
            # Thêm thư mục proto của bạn vào đường dẫn tìm kiếm
            f"-I{PROTO_DIR}",
            # File proto nguồn cần xử lý
            proto_file_path
        ]
        # --- KẾT THÚC PHẦN SỬA LỖI ---
        
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            print(f"Successfully generated files for {filename}")
        except subprocess.CalledProcessError as e:
            print(f"Error generating files for {filename}:")
            # In ra lỗi từ stderr để gỡ lỗi dễ hơn
            print(e.stderr)
            env.Exit(1)

print("Protobuf generation finished.")