# ===============================================================================
# --------------------------- INFORMATIONS GENERALES ----------------------------
# ===============================================================================
# C++ standard
set(CMAKE_C_COMPILER    "/opt/gcc-8.3.0/bin/gcc")
set(CMAKE_CXX_COMPILER  "/opt/gcc-8.3.0/bin/g++")
set (CMAKE_CXX_STANDARD "17")

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# Définition du projet et de sa version
project(rk_games_project VERSION 1.0)
set(EXE_NAME "rk_games")

# Fichier de configuration interne
include(GNUInstallDirs)
configure_file(src/config.h.in config.h)

# ===============================================================================
# ------------------------------ FICHIERS SOURCES -------------------------------
# ===============================================================================
FILE(GLOB_RECURSE SOURCE_FILES src/*.cpp include/*cpp)

# ===============================================================================
# ----------------------------- INCLUDE DIRECTORIES -----------------------------
# ===============================================================================
include_directories(include)
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE
                           ${PROJECT_BINARY_DIR}    
                           ${PROJECT_SOURCE_DIR}/src>
                           )

# ===============================================================================
# -------------------------------- LIBRAIRIES -----------------------------------
# ===============================================================================

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---- Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
     target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
endif()

# ---- OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# ---- Boost
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS iostreams system filesystem) 
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})  
    target_link_libraries(${PROJECT_NAME} INTERFACE ${Boost_LIBRARIES})
endif()

# ----- Matplotlib
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE
  Python3::Python
)
find_package(Python3 COMPONENTS NumPy)
if(Python3_NumPy_FOUND)
  target_link_libraries(${PROJECT_NAME} INTERFACE
    Python3::NumPy
  )
else()
  target_compile_definitions(${PROJECT_NAME} INTERFACE WITHOUT_NUMPY)
endif()

# ===============================================================================
# -------------------------- CREATION DE L'EXECUTABLE ---------------------------
# ===============================================================================

add_executable(${EXE_NAME} ${SOURCE_FILES})
target_link_libraries(${EXE_NAME} PRIVATE ${PROJECT_NAME})

# ===============================================================================
# ---------------------------- OPTIONS DE COMPILATION ---------------------------
# ===============================================================================

#set(CMAKE_DEBUG_POSTFIX "_debug")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/build")
#set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} /MD /O2 -O3")

target_compile_options(${EXE_NAME} PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Werror -O3 -Wno-register>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W0 /O2>)
