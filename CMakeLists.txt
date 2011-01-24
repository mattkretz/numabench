cmake_minimum_required(VERSION 2.6.2)
project(Vc)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$")
   set(COMPILER_IS_INTEL true)
   set(COMPILER_IS_OPEN64 false)
   set(COMPILER_IS_MSVC false)
   set(COMPILER_IS_GCC false)
elseif(CMAKE_CXX_COMPILER MATCHES "/(opencc|openCC)$")
   set(COMPILER_IS_INTEL false)
   set(COMPILER_IS_OPEN64 true)
   set(COMPILER_IS_MSVC false)
   set(COMPILER_IS_GCC false)
elseif(MSVC)
   set(COMPILER_IS_INTEL false)
   set(COMPILER_IS_OPEN64 false)
   set(COMPILER_IS_MSVC true)
   set(COMPILER_IS_GCC false)
elseif(CMAKE_COMPILER_IS_GNUCXX)
   set(COMPILER_IS_INTEL false)
   set(COMPILER_IS_OPEN64 false)
   set(COMPILER_IS_MSVC false)
   set(COMPILER_IS_GCC true)
else()
   message(ERROR "Unsupported Compiler. Please fill out the missing parts in the CMake scripts and submit a patch to http://code.compeng.uni-frankfurt.de/projects/vc")
endif()

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include (MacroEnsureVersion)
include (OptimizeForArchitecture)
include (AddCompilerFlag)
include (CheckCXXSourceRuns)
include (AddTargetProperty)

check_cxx_source_runs("int main() { return sizeof(void*) != 8; }" VOID_PTR_IS_64BIT)

if(NO_BSF_LOOPS)
   add_definitions(-DVC_NO_BSF_LOOPS)
endif(NO_BSF_LOOPS)
if(NO_GATHER_TRICKS)
   add_definitions(-DVC_NO_GATHER_TRICKS)
endif(NO_GATHER_TRICKS)
if(GATHER_SET)
   add_definitions(-DVC_GATHER_SET)
endif(GATHER_SET)

set(SSE_INTRINSICS_BROKEN false)

if(COMPILER_IS_OPEN64)
   # Open64 is detected as GNUCXX :(
   AddCompilerFlag("-W")
   AddCompilerFlag("-Wall")
   AddCompilerFlag("-Wimplicit")
   AddCompilerFlag("-Wswitch")
   AddCompilerFlag("-Wformat")
   AddCompilerFlag("-Wchar-subscripts")
   AddCompilerFlag("-Wparentheses")
   AddCompilerFlag("-Wmultichar")
   AddCompilerFlag("-Wtrigraphs")
   AddCompilerFlag("-Wpointer-arith")
   AddCompilerFlag("-Wcast-align")
   AddCompilerFlag("-Wreturn-type")
   AddCompilerFlag("-Wno-unused-function")
   AddCompilerFlag("-ansi")
   AddCompilerFlag("-pedantic")
   AddCompilerFlag("-Wno-long-long")
   AddCompilerFlag("-Wshadow")
   AddCompilerFlag("-Wold-style-cast")
   AddCompilerFlag("-Wno-variadic-macros")
   AddCompilerFlag("-fno-threadsafe-statics")
   set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
   set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ")
   set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
   set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
   set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL} ")
   set(CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE} -O3")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")

   # if compiling for 32 bit x86 we need to use the -mfpmath=sse since the x87 is broken by design
   if(NOT VOID_PTR_IS_64BIT)
      exec_program(${CMAKE_C_COMPILER} ARGS -dumpmachine OUTPUT_VARIABLE _gcc_machine)
      if(_gcc_machine MATCHES "[x34567]86")
         AddCompilerFlag("-mfpmath=sse")
      endif(_gcc_machine MATCHES "[x34567]86")
   endif(NOT VOID_PTR_IS_64BIT)

   # Open64 uses binutils' assembler. And it has to be recent enough otherwise it'll bail out at some of the instructions Vc uses
   find_program(_binutils_as as)
   if(NOT _binutils_as)
      message(WARNING "Could not find 'as', the assembly normally used by GCC. Hoping everything will work out...")
   else(NOT _binutils_as)
      if(APPLE)
         # it's not really binutils, but it'll give us the assembler version which is what we want
         exec_program(${_binutils_as} ARGS -v /dev/null OUTPUT_VARIABLE _as_version)
      else(APPLE)
         exec_program(${_binutils_as} ARGS --version OUTPUT_VARIABLE _as_version)
      endif(APPLE)
      string(REGEX REPLACE "\\([^\\)]*\\)" "" _as_version "${_as_version}")
      string(REGEX MATCH "[1-9]\\.[0-9]+(\\.[0-9]+)?" _as_version "${_as_version}")
      macro_ensure_version("2.18.93" "${_as_version}" _as_good)
      if(NOT _as_good)
         message(SEND_ERROR "Your binutils is too old (${_as_version}). The assembler will not be able to compile Vc.")
      endif(NOT _as_good)
   endif(NOT _binutils_as)
   mark_as_advanced(_binutils_as)
elseif(COMPILER_IS_GCC)
   ##################################################################################################
   #                                              GCC                                               #
   ##################################################################################################
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -W -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function -ansi -pedantic -Wno-long-long -Wshadow")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function -ansi -pedantic -Wno-long-long -Wshadow")
   AddCompilerFlag("-Wold-style-cast")
   AddCompilerFlag("-Wno-variadic-macros")
   AddCompilerFlag("-fno-threadsafe-statics")
   set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}")
   set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} ")
   set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
   set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
   set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL} ")
   set(CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE} -O3")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")

   # check the GCC version
   exec_program(${CMAKE_C_COMPILER} ARGS -dumpversion OUTPUT_VARIABLE _gcc_version)
   macro_ensure_version("4.4.1" "${_gcc_version}" GCC_4_4_1)
   if(NOT GCC_4_4_1)
      macro_ensure_version("4.3.0" "${_gcc_version}" GCC_4_3_0)
      if(NOT GCC_4_3_0)
         message(STATUS "WARNING: Your GCC is older than 4.3.0. It is unable to handle all SSE2 intrinsics. All SSE code will be disabled. Please update to the latest GCC if you can.\n-- \n-- ")
         set(SSE_INTRINSICS_BROKEN true)
      endif(NOT GCC_4_3_0)
   endif(NOT GCC_4_4_1)

   # if compiling for 32 bit x86 we need to use the -mfpmath=sse since the x87 is broken by design
   if(NOT VOID_PTR_IS_64BIT)
      exec_program(${CMAKE_C_COMPILER} ARGS -dumpmachine OUTPUT_VARIABLE _gcc_machine)
      if(_gcc_machine MATCHES "[x34567]86")
         AddCompilerFlag("-mfpmath=sse")
      endif(_gcc_machine MATCHES "[x34567]86")
   endif(NOT VOID_PTR_IS_64BIT)

   # GCC uses binutils' assembler. And it has to be recent enough otherwise it'll bail out at some of the instructions Vc uses
   find_program(_binutils_as as)
   if(NOT _binutils_as)
      message(WARNING "Could not find 'as', the assembler normally used by GCC. Hoping everything will work out...")
   else(NOT _binutils_as)
      if(APPLE)
         # it's not really binutils, but it'll give us the assembler version which is what we want
         exec_program(${_binutils_as} ARGS -v /dev/null OUTPUT_VARIABLE _as_version)
      else(APPLE)
         exec_program(${_binutils_as} ARGS --version OUTPUT_VARIABLE _as_version)
      endif(APPLE)
      string(REGEX REPLACE "\\([^\\)]*\\)" "" _as_version "${_as_version}")
      string(REGEX MATCH "[1-9]\\.[0-9]+(\\.[0-9]+)?" _as_version "${_as_version}")
      macro_ensure_version("2.18.93" "${_as_version}" _as_good)
      if(NOT _as_good)
         message(WARNING "Your binutils is too old (${_as_version}). Some optimizations of Vc will be disabled.")
         add_definitions(-DVC_NO_XGETBV) # old assembler doesn't know the xgetbv instruction
      endif(NOT _as_good)
   endif(NOT _binutils_as)
   mark_as_advanced(_binutils_as)
elseif(COMPILER_IS_INTEL)
   ##################################################################################################
   #                                          Intel Compiler                                        #
   ##################################################################################################

   set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
   set(CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE} -O3")
   set(CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")

   set(ALIAS_FLAGS "-no-ansi-alias")

   # iomanip from latest libstdc++ makes ICC fail unless C++0x is selected
   AddCompilerFlag("-std=c++0x")

   # per default icc is not IEEE compliant, but we need that for verification
   set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${ALIAS_FLAGS} -w1 -fp-model precise")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ALIAS_FLAGS} -w1 -fp-model precise")
elseif(COMPILER_IS_MSVC)
   AddCompilerFlag("/wd4800") # Disable warning "forcing value to bool"
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)
   # MSVC does not support inline assembly on 64 bit! :(
   # searching the help for xgetbv doesn't turn up anything. So just fall back to not supporting AVX on Windows :(
   add_definitions(-DVC_NO_XGETBV)
endif()

OptimizeForArchitecture()

if(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
   message(STATUS "WARNING! It seems you are compiling without optimization. Please set CMAKE_BUILD_TYPE.")
endif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Vc STATIC vector.cpp common/support.cpp)
add_target_property(Vc COMPILE_FLAGS "-DVC_COMPILE_LIB")

find_file(HAVE_SYS_MMAN sys/mman.h)
if(NOT HAVE_SYS_MMAN)
   message(ERROR "could not find required sys/mman.h header to enable memlocking")
   add_definitions(-DVC_BENCHMARK_NO_MLOCK)
endif()
mark_as_advanced(HAVE_SYS_MMAN)

add_library(cpuset STATIC benchmarks/cpuset.c)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   set(REALTIME_BENCHMARKS FALSE CACHE BOOL "Use SCHED_FIFO priority to give benchmarks highest priority, achieving less jitter in the results.")
   if(REALTIME_BENCHMARKS)
      add_executable(benchmark benchmarks/benchmark.cpp)
      add_definitions(-DSCHED_FIFO_BENCHMARKS)
      add_custom_command(TARGET benchmark POST_BUILD
         COMMAND sudo chown root benchmark
         COMMAND sudo chmod +s benchmark
         WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
         COMMENT "Making the benchmark helper suid root"
         )
   endif(REALTIME_BENCHMARKS)
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

add_executable(numabench benchmarks/numabench.cpp)
add_target_property(numabench COMPILE_FLAGS "-DVC_IMPL=SSE")
target_link_libraries(numabench Vc cpuset)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   target_link_libraries(numabench rt)
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

configure_file(processOutput.sh processOutput.sh COPYONLY)
