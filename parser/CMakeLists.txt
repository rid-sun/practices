cmake_minimum_required(VERSION 3.0.0)
project(parser VERSION 0.1.0)

# 创建out并复制输出
set(out ${CMAKE_CURRENT_SOURCE_DIR}/bin)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${out})

# 注册子级别目录
add_subdirectory(src)