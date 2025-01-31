cmake_minimum_required (VERSION 2.6)
project (TranslateAllLOLA)

set (CMAKE_CXX_STANDARD 11)

# add the binary tree to the search path for include files
# so that we will find TestSerialConfig.h
include_directories ("${PROJECT_BINARY_DIR}")

# add the executable
add_executable (TranslateAllLOLA TranslateAllLOLA.cpp)
target_link_libraries (TranslateAllLOLA  ${EXTRA_LIBS})

# add the install targets
install (TARGETS TranslateAllLOLA DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIR}/TranslateAllLOLA.h"
  DESTINATION include)
