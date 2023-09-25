add_library(nsga STATIC ./nsga.cpp ./nsga.hpp)
target_include_directories(nsga PUBLIC . )
set_property(TARGET nsga PROPERTY CXX_STANDARD 20)
# 设置公共头文件的属性
set_target_properties(nsga
    PROPERTIES
    PUBLIC_HEADER ./nsga.hpp # 这里再次指定公共头文件
)