cmake_minimum_required(VERSION 3.10)
project(senkeidaisu)

# Устанавливаем стандарт C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Устанавливаем пути к исходникам и заголовочным файлам
include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/include)

# Добавляем все исходные файлы из src/
file(GLOB SOURCES src/*.hpp src/*.cpp)

# Создание библиотеки
add_library(senkeidaisu_lib ${SOURCES})

# Исполняемый файл (если нужен)
add_executable(senkeidaisu src/main.cpp)
target_link_libraries(senkeidaisu senkeidaisu_lib)

# Добавляем тесты (если будут)
# add_subdirectory(tests)
