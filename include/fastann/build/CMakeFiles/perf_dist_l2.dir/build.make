# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ncarrara/workspace/EDMC/include/fastann

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ncarrara/workspace/EDMC/include/fastann/build

# Include any dependencies generated for this target.
include CMakeFiles/perf_dist_l2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/perf_dist_l2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/perf_dist_l2.dir/flags.make

CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.o: CMakeFiles/perf_dist_l2.dir/flags.make
CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.o: ../perf_dist_l2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ncarrara/workspace/EDMC/include/fastann/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.o -c /home/ncarrara/workspace/EDMC/include/fastann/perf_dist_l2.cpp

CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ncarrara/workspace/EDMC/include/fastann/perf_dist_l2.cpp > CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.i

CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ncarrara/workspace/EDMC/include/fastann/perf_dist_l2.cpp -o CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.s

# Object files for target perf_dist_l2
perf_dist_l2_OBJECTS = \
"CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.o"

# External object files for target perf_dist_l2
perf_dist_l2_EXTERNAL_OBJECTS =

perf_dist_l2: CMakeFiles/perf_dist_l2.dir/perf_dist_l2.cpp.o
perf_dist_l2: CMakeFiles/perf_dist_l2.dir/build.make
perf_dist_l2: libfastann.so
perf_dist_l2: CMakeFiles/perf_dist_l2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ncarrara/workspace/EDMC/include/fastann/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable perf_dist_l2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/perf_dist_l2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/perf_dist_l2.dir/build: perf_dist_l2

.PHONY : CMakeFiles/perf_dist_l2.dir/build

CMakeFiles/perf_dist_l2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/perf_dist_l2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/perf_dist_l2.dir/clean

CMakeFiles/perf_dist_l2.dir/depend:
	cd /home/ncarrara/workspace/EDMC/include/fastann/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ncarrara/workspace/EDMC/include/fastann /home/ncarrara/workspace/EDMC/include/fastann /home/ncarrara/workspace/EDMC/include/fastann/build /home/ncarrara/workspace/EDMC/include/fastann/build /home/ncarrara/workspace/EDMC/include/fastann/build/CMakeFiles/perf_dist_l2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/perf_dist_l2.dir/depend

