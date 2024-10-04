# 使用 LKH2, 要用 LKH3 的话就把下面所有出现的 LKH2 改成 LKH3

# cmake_minimum_required(VERSION 2.8)
cmake_minimum_required(VERSION 3.10.2)
project(project_hw)

set(CMAKE_C_STANDARD 99)

# 设置路径变量
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(ALG_DIR ${PROJECT_SOURCE_DIR}/algorithm)
set(LKH_DIR ${PROJECT_SOURCE_DIR}/LKH3)

# 设置头文件目录变量
set(HEADER_DIRS
    ${ALG_DIR}
    ${LKH_DIR}
    ${LKH_DIR}/include
)

# 添加源文件变量
set(SOURCES
    ${PROJECT_SOURCE_DIR}/main.c
    ${ALG_DIR}/algorithm.c
    ${LKH_DIR}/LKHInterface.c
)

# 包含头文件目录
include_directories(${HEADER_DIRS})

# LKH 库选项
# ------------------ option 1: from src ---------------------
# 从 LKH 源码 编译成 LKH2 库 或 LKH3 库
aux_source_directory(${LKH_DIR}/src LKH_SRC)
add_library(LKH3 ${LKH_SRC})
# ------------------ option 1: from src ---------------------

# --------------------- option 2: from lib ---------------------
# 直接链接之前编译好的 LKH2 库或 LKH3 库，其路径由下文 link_directories指定
# --------------------- option 2: from lib ---------------------

# 设置可执行文件输出路径
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# 指定库文件夹
link_directories(${PROJECT_SOURCE_DIR}/lib)

# 生成可执行文件
add_executable(project_hw  ${SOURCES})

# 链接库 libseek_model.so libLKH2.a libm.so(数学库)
target_link_libraries(project_hw PRIVATE seek_model LKH3 m)

# 添加编译选项和宏定义
target_compile_options(project_hw PRIVATE -g -Wall -std=c99)
target_compile_definitions(project_hw PRIVATE _GNU_SOURCE)
