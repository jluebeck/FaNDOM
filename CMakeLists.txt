cmake_minimum_required(VERSION 3.1)

#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package (Threads REQUIRED)

project(FaNDOM)

set(CMAKE_CXX_STANDARD 11)

add_executable(FaNDOM main.cpp OMIO.cpp OMIO.h OMHelper.cpp OMHelper.h threadsafe_queue.h fandomAligner.cpp fandomAligner.h OMFilter.cpp OMFilter.h SeedGraph.cpp SeedGraph.h)

#target_link_libraries(FaNDOM Threads::Threads)
TARGET_LINK_LIBRARIES(FaNDOM pthread)
