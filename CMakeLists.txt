cmake_minimum_required(VERSION 3.13)
project(vslam)

set(CMAKE_CXX_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)



find_package(OpenCV 3 REQUIRED)
find_package(Eigen3)
find_package(Sophus REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${OpenCV_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${Sophus_INCLUDE_DIRS}
)

set(THIRD_PARTY_LIBS
        ${OpenCV_LIBS}
        ${EIGEN3_LIBS}
        ${Sophus_LIBRARIES}
)

add_library(${PROJECT_NAME} SHARED
        src/camera.cpp
        src/config.cpp
#        src/mappoint.cpp
#        src/frame.cpp
#        src/map.cpp
        src/bundleadjustment.cpp)

target_link_libraries(${PROJECT_NAME} ${THIRD_PARTY_LIBS})

add_executable(run_slam test/run_slam.cpp)
#add_executable(ttt test/ttt.cpp)

target_link_libraries(run_slam ${PROJECT_NAME})