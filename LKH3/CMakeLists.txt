# cmake_minimum_required(VERSION 2.8)
cmake_minimum_required(VERSION 3.10.2)
Project(testLKH3_Interface)

set(CMAKE_C_STANDARD 99)

# 设置路径变量
set(LKH_SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(LKH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# 设置头文件目录变量
set(HEADER_DIRS
    ${LKH_INCLUDE_DIR}
)

# 添加源文件变量
set(SOURCES
    ${PROJECT_SOURCE_DIR}/testLKH.c
    ${PROJECT_SOURCE_DIR}/LKHInterface.c
)

# 包含头文件目录
include_directories(${HEADER_DIRS})

# 编译 LKH 3 源码成 LKH3 库
aux_source_directory(${LKH_SRC_DIR} LKH_SRC)
add_library(LKH3 ${LKH_SRC})

# 设置可执行文件输出路径
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# 生成可执行文件
add_executable(testLKH ${SOURCES})

# 链接库
target_link_libraries(testLKH LKH3 m) # 链接 LKH3 库和数学库
 
# 添加编译选项和宏定义
target_compile_options(testLKH PRIVATE -g -Wall -std=c99)
target_compile_definitions(testLKH PRIVATE _GNU_SOURCE)
