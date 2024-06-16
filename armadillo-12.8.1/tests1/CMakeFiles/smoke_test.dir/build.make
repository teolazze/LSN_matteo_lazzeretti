# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1

# Include any dependencies generated for this target.
include tests1/CMakeFiles/smoke_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests1/CMakeFiles/smoke_test.dir/compiler_depend.make

# Include the progress variables for this target.
include tests1/CMakeFiles/smoke_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests1/CMakeFiles/smoke_test.dir/flags.make

tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.o: tests1/CMakeFiles/smoke_test.dir/flags.make
tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.o: tests1/smoke_test.cpp
tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.o: tests1/CMakeFiles/smoke_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.o"
	cd /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.o -MF CMakeFiles/smoke_test.dir/smoke_test.cpp.o.d -o CMakeFiles/smoke_test.dir/smoke_test.cpp.o -c /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1/smoke_test.cpp

tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smoke_test.dir/smoke_test.cpp.i"
	cd /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1/smoke_test.cpp > CMakeFiles/smoke_test.dir/smoke_test.cpp.i

tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smoke_test.dir/smoke_test.cpp.s"
	cd /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1/smoke_test.cpp -o CMakeFiles/smoke_test.dir/smoke_test.cpp.s

# Object files for target smoke_test
smoke_test_OBJECTS = \
"CMakeFiles/smoke_test.dir/smoke_test.cpp.o"

# External object files for target smoke_test
smoke_test_EXTERNAL_OBJECTS =

tests1/smoke_test: tests1/CMakeFiles/smoke_test.dir/smoke_test.cpp.o
tests1/smoke_test: tests1/CMakeFiles/smoke_test.dir/build.make
tests1/smoke_test: libarmadillo.so.12.8.1
tests1/smoke_test: /usr/lib/x86_64-linux-gnu/libblas.so
tests1/smoke_test: /usr/lib/x86_64-linux-gnu/liblapack.so
tests1/smoke_test: /usr/lib/x86_64-linux-gnu/libarpack.so
tests1/smoke_test: /usr/lib/x86_64-linux-gnu/libsuperlu.so
tests1/smoke_test: tests1/CMakeFiles/smoke_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable smoke_test"
	cd /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/smoke_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests1/CMakeFiles/smoke_test.dir/build: tests1/smoke_test
.PHONY : tests1/CMakeFiles/smoke_test.dir/build

tests1/CMakeFiles/smoke_test.dir/clean:
	cd /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 && $(CMAKE_COMMAND) -P CMakeFiles/smoke_test.dir/cmake_clean.cmake
.PHONY : tests1/CMakeFiles/smoke_test.dir/clean

tests1/CMakeFiles/smoke_test.dir/depend:
	cd /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1 /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1 /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1 /home/matteo/Documenti/LSN2024/NSL_SIMULATOR/armadillo-12.8.1/tests1/CMakeFiles/smoke_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests1/CMakeFiles/smoke_test.dir/depend

