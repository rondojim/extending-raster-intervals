cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)

# Add Boost
find_package(Boost REQUIRED COMPONENTS)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})

# Include directories
include_directories(include)

# Source files (excluding main.cc and produce_hilbert.cc)
file(GLOB_RECURSE SOURCE_FILES
    "src/utils/*.cc"
    "src/join_algos/*.cc"
    "src/mbr_algos/*.cc"
)

# Main executable
add_executable(main_executable src/main.cc ${SOURCE_FILES})
target_link_libraries(main_executable ${Boost_LIBRARIES})

# Hilbert executable
add_executable(hilbert_executable src/produce_hilbert.cc ${SOURCE_FILES})
target_link_libraries(hilbert_executable ${Boost_LIBRARIES})

# Tests
set(TEST_FILES
    tests/join_algos/test_ri_join_op.cc
    tests/mbr_algos/test_forward_scan.cc
    # tests/join_algos/test_mbr_filter.cc
    # tests/utils/test_data_reader.cc
    # tests/utils/test_geometry_types.cc
    # tests/utils/test_hilbert.cc
)

# Add test executables
foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE} ${SOURCE_FILES})
    target_link_libraries(${TEST_NAME} ${Boost_LIBRARIES})
endforeach()

# Enable testing
enable_testing()

# Add tests
foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
