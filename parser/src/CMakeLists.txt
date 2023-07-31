# 注意，在这里键入你的python安装包的include文件夹路径
include_directories("D:/Anaconda3/include")
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/eigen-3.4.0)

# 注意，在这里键入你的python安装包的libs文件夹路径
link_directories("D:/Anaconda3/libs") 

# 寻找所有的src文件夹下的cpp文件
file(GLOB srcs ./*.cpp)

# 生成可执行文件parser
add_executable(parser ${srcs})
target_link_libraries(parser python39.lib) # 注意，在这里要链接静态库

# 指定parser输出路径
set_target_properties(parser PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
