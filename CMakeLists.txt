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
add_library(EFDLib STATIC EFD.cpp)
target_link_libraries(EFDLib PUBLIC OptionType)  # PUBLIC allows EFD dependents to access UtilLib

# BSM Library
add_library(BSMLib STATIC BSM.cpp)
target_link_libraries(BSMLib PUBLIC UtilLib)

# IFD Library
add_library(IFDLib STATIC IFD.cpp)
target_link_libraries(IFDLib PUBLIC OptionType)

# CNFD Library
add_library(CNFDLib STATIC CNFD.cpp)
target_link_libraries(CNFDLib PUBLIC OptionType)

# Option Data Library
add_library(OptionDataLib STATIC OptionData.cpp)
target_link_libraries(OptionDataLib PUBLIC BSMLib EFDLib IFDLib CNFDLib)

# Ticker Library
add_library(TickerLib STATIC Ticker.cpp)
target_link_libraries(TickerLib PUBLIC OptionDataLib)

# FunctionsUtil Library
add_library(FunctionsUtilLib STATIC FunctionsUtil.cpp)
target_link_libraries(FunctionsUtilLib PUBLIC TickerLib)

# Main executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC FunctionsUtilLib)
