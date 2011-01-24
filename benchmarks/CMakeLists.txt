find_file(HAVE_SYS_MMAN sys/mman.h)
if(HAVE_SYS_MMAN)
   set(USE_MLOCK_IN_BENCHMARKS TRUE CACHE BOOL "Use memlock to get more reliable results in the benchmarks.")
   if(NOT USE_MLOCK_IN_BENCHMARKS)
      add_definitions(-DVC_BENCHMARK_NO_MLOCK)
   endif(NOT USE_MLOCK_IN_BENCHMARKS)
else(HAVE_SYS_MMAN)
   add_definitions(-DVC_BENCHMARK_NO_MLOCK)
endif(HAVE_SYS_MMAN)
mark_as_advanced(HAVE_SYS_MMAN)

add_library(cpuset STATIC cpuset.c)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   set(REALTIME_BENCHMARKS FALSE CACHE BOOL "Use SCHED_FIFO priority to give benchmarks highest priority, achieving less jitter in the results.")
   if(REALTIME_BENCHMARKS)
      add_executable(benchmark benchmark.cpp)
      add_definitions(-DSCHED_FIFO_BENCHMARKS)
      add_custom_command(TARGET benchmark POST_BUILD
         COMMAND sudo chown root benchmark
         COMMAND sudo chmod +s benchmark
         WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
         COMMENT "Making the benchmark helper suid root"
         )
   endif(REALTIME_BENCHMARKS)
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

add_executable(membench_sse membench.cpp)
add_target_property(membench_sse COMPILE_FLAGS "-DVC_IMPL=SSE")
target_link_libraries(membench_sse Vc cpuset)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   target_link_libraries(membench_sse rt)
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
