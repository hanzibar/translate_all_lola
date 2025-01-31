cmake_minimum_required (VERSION 3.10)
project (TranslateAllLOLA)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# add the binary tree to the search path for include files
# so that we will find TestSerialConfig.h
include_directories ("${PROJECT_BINARY_DIR}")

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})

# add the executable
add_executable (TranslateAllLOLA 
    TranslateAllLOLA.cpp
    libbmp.cpp
)

target_link_libraries (TranslateAllLOLA ${PNG_LIBRARIES})

# add the install targets
install (TARGETS TranslateAllLOLA DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIR}/TranslateAllLOLA.h"
  DESTINATION include)
