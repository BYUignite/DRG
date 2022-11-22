cmake_minimum_required(VERSION 3.1)

project(drg LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(drg.x)
target_sources(drg.x PRIVATE driver.cc drg.cc streams.cc)

find_library(KINSOL sundials_kinsol)
find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)
find_package(Python COMPONENTS Interpreter Development)
if(APPLE)
    find_library(ACCELERATE Accelerate)
    target_include_directories(drg.x PRIVATE "/opt/homebrew/include/")
endif()

set(CMAKE_BUILD_RPATH "/usr/local/lib")

target_include_directories(drg.x PRIVATE "/usr/local/include")
target_link_libraries(drg.x ${CANTERA} ${KINSOL} ${CVODES} ${Python_LIBRARIES} ${ACCELERATE})

install(TARGETS drg.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
