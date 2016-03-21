cmake_minimum_required(VERSION 3.2)
project(NetworkFlow)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.54.0 COMPONENTS program_options)

set(gtest_SOURCE_DIR /media/alvis/Data/MasterThesis/NetworkFlow/lib/googletest/googletest)

include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES Graph.cpp Graph.h Common.h parallel.h SIA.cpp SIA.h TimerTool.cpp TimerTool.h
        CostScaling.cpp CostScaling.h NodeList.h LocalDominant.cpp LocalDominant.h Lemon.h SCS.cpp SCS.h)
add_executable(NetworkFlows main.cpp ${SOURCE_FILES})
add_executable(Testing testing.cpp ${SOURCE_FILES})
add_executable(Generator generator.cpp ${SOURCE_FILES})

target_link_libraries(NetworkFlows ${Boost_LIBRARIES})
target_link_libraries(Generator ${Boost_LIBRARIES})
target_link_libraries(Testing ${Boost_LIBRARIES})