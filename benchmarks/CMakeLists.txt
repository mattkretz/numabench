include(BenchmarkMacros)

set(NO_AUTOVEC "-fno-tree-vectorize")
check_cxx_compiler_flag(${NO_AUTOVEC} check_compiler_flag_no_autovect_gcc)
if(NOT check_compiler_flag_no_autovect_gcc)
   set(NO_AUTOVEC "-no-vec")
   check_cxx_compiler_flag(${NO_AUTOVEC} check_compiler_flag_no_autovect_icc)
   if(NOT check_compiler_flag_no_autovect_icc)
      message(STATUS "Missing flag to disable autovectorization in the compiler")
      SET(NO_AUTOVEC)
   endif(NOT check_compiler_flag_no_autovect_icc)
endif(NOT check_compiler_flag_no_autovect_gcc)

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

vc_add_benchmark(memio)
add_target_property(memio_scalar_benchmark COMPILE_FLAGS ${NO_AUTOVEC})

add_executable(membench_sse membench.cpp)
add_target_property(membench_sse COMPILE_FLAGS "-DVC_IMPL=SSE")
target_link_libraries(membench_sse Vc ${LIBS})
