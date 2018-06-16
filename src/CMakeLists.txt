add_subdirectory(test/)
add_subdirectory(util/)

add_executable(partition partition.cpp)
target_link_libraries(partition util ${CMAKE_THREAD_LIBS_INIT} "${CMAKE_SOURCE_DIR}/processorcountermonitor/libPCM.a" papi numa)
