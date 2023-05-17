cmake_minimum_required(VERSION 3.0.0)
project(parser VERSION 0.1.0)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories("D:/Anaconda3/include") # 注意，在这里键入你的python安装包的include文件夹路径
include_directories(include)

# 注意，在这里键入你的python安装包的libs文件夹路径
link_directories("D:/Anaconda3/libs") 

file(GLOB srcs src/*.cpp)
add_executable(parser main.cpp ${srcs})
target_link_libraries(parser python39.lib) # 注意，在这里要链接静态库

# 创建out并复制输出
set(out ${CMAKE_CURRENT_SOURCE_DIR}/out)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${out})
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/build/parser.exe ${out})