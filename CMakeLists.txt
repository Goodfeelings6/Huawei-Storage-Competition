# cmake_minimum_required(VERSION 2.8)
cmake_minimum_required(VERSION 3.10.2)
project(project_hw)

set(CMAKE_C_STANDARD 99)

# 设置全局编译选项
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -std=c99")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -std=c99")

# 或者指定构建方式为 release 模式
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

# 设置路径变量
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(ALG_DIR ${PROJECT_SOURCE_DIR}/algorithm)
set(LKH_DIR ${PROJECT_SOURCE_DIR}/LKH2_MIN)

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
# 从 LKH 源码 编译成 LKH 库
aux_source_directory(${LKH_DIR}/src LKH_SRC)
add_library(LKH ${LKH_SRC})
# ------------------ option 1: from src ---------------------

# --------------------- option 2: from lib ---------------------
# 直接链接之前编译好的 LKH 库，其路径由下文 link_directories指定
# --------------------- option 2: from lib ---------------------

# 设置可执行文件输出路径
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# 指定库文件夹
link_directories(${PROJECT_SOURCE_DIR}/lib)

# 生成可执行文件
add_executable(project_hw  ${SOURCES})

# 链接库 libseek_model.so libLKH.a libm.so(数学库)
target_link_libraries(project_hw PRIVATE seek_model LKH m)

# 添加编译选项和宏定义
# target_compile_options(project_hw PRIVATE -g -O2 -Wall -std=c99)
target_compile_definitions(project_hw PRIVATE _GNU_SOURCE)
