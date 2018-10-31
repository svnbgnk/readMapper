# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (readMapper CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
#find_package (ZLIB)

set (CMAKE_PREFIX_PATH "/home/sven/devel/seqan3/build_system/")
#set (SEQAN_INCLUDE_PATH "/home/sven/devel/seqan3/include")

# Load the SeqAn module and fail if not found.
find_package (SeqAn3 REQUIRED)
find_package (OpenMP)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN3_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN3_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS " -std=c++17 -fconcepts -DSEQAN_DISABLE_VERSION_CHECK=YES ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (readMapper readMapper.cpp)
target_link_libraries (readMapper ${SEQAN3_LIBRARIES})

add_executable (create_index create_index.cpp)
target_link_libraries (create_index ${SEQAN3_LIBRARIES})

add_executable (mappability mappability.cpp)
target_link_libraries (mappability ${SEQAN3_LIBRARIES})

add_executable (create_bit_v create_bit_v.cpp)
target_link_libraries (create_bit_v ${SEQAN3_LIBRARIES})
