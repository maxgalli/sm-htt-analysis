set(current_dir ${CMAKE_SOURCE_DIR}/jsoncpp)

add_library(jsoncpp SHARED ${current_dir}/jsoncpp.cpp)
target_include_directories(jsoncpp PRIVATE ${current_dir}/json)

