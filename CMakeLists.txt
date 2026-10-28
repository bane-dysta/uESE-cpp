cmake_minimum_required(VERSION 3.16)
project(uese_cpp VERSION 1.1 LANGUAGES CXX)

add_library(uese_core STATIC
    src/io.cpp src/model.cpp src/output.cpp src/parameters.cpp src/surface.cpp)
target_include_directories(uese_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(uese_core PUBLIC cxx_std_17)
set_target_properties(uese_core PROPERTIES CXX_EXTENSIONS OFF)

add_executable(uese src/main.cpp)
target_link_libraries(uese PRIVATE uese_core)
if(MSVC)
    target_compile_options(uese_core PRIVATE /W4)
else()
    target_compile_options(uese_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME uese_regression COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py $<TARGET_FILE:uese>)
    endif()
endif()

install(TARGETS uese RUNTIME DESTINATION bin)
