cmake_minimum_required(VERSION 3.16.0)

project(test VERSION 0.1.0 LANGUAGES C CXX)

add_subdirectory(src/nsga)
add_executable(example ${PROJECT_SOURCE_DIR}/example/main.cpp)
target_link_libraries(example PRIVATE nsga)

