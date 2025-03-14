cmake_minimum_required(VERSION 3.16)
project(FE621HW2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library for OptionType
add_library(OptionType INTERFACE)
target_include_directories(OptionType INTERFACE ${CMAKE_SOURCE_DIR})

# Utility Library (Common Functions)
add_library(UtilLib STATIC util.cpp)
target_link_libraries(UtilLib PUBLIC OptionType)  # PUBLIC allows all users of UtilLib to access OptionType

# EFD Library
add_library(EFDlib STATIC EFD.cpp)
target_link_libraries(EFDlib PUBLIC OptionType)  # PUBLIC allows EFD dependents to access UtilLib

# BSM Library
add_library(BSMlib STATIC BSM.cpp)
target_link_libraries(BSMlib PUBLIC UtilLib)


# Main executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC EFDlib BSMlib)
