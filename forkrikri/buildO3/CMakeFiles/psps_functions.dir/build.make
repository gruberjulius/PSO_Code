# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cluster/home/pafloria/testMoe5/Moe

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cluster/home/pafloria/testMoe5/buildO3

# Include any dependencies generated for this target.
include CMakeFiles/psps_functions.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/psps_functions.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/psps_functions.dir/flags.make

CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o: CMakeFiles/psps_functions.dir/flags.make
CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o: /cluster/home/pafloria/testMoe5/Moe/examples/psps_functions.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /cluster/home/pafloria/testMoe5/buildO3/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o -c /cluster/home/pafloria/testMoe5/Moe/examples/psps_functions.cpp

CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -E /cluster/home/pafloria/testMoe5/Moe/examples/psps_functions.cpp > CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.i

CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -S /cluster/home/pafloria/testMoe5/Moe/examples/psps_functions.cpp -o CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.s

CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.requires:
.PHONY : CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.requires

CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.provides: CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.requires
	$(MAKE) -f CMakeFiles/psps_functions.dir/build.make CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.provides.build
.PHONY : CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.provides

CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.provides.build: CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o

# Object files for target psps_functions
psps_functions_OBJECTS = \
"CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o"

# External object files for target psps_functions
psps_functions_EXTERNAL_OBJECTS =

psps_functions: CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o
psps_functions: CMakeFiles/psps_functions.dir/build.make
psps_functions: CMakeFiles/psps_functions.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable psps_functions"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/psps_functions.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/psps_functions.dir/build: psps_functions
.PHONY : CMakeFiles/psps_functions.dir/build

CMakeFiles/psps_functions.dir/requires: CMakeFiles/psps_functions.dir/examples/psps_functions.cpp.o.requires
.PHONY : CMakeFiles/psps_functions.dir/requires

CMakeFiles/psps_functions.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/psps_functions.dir/cmake_clean.cmake
.PHONY : CMakeFiles/psps_functions.dir/clean

CMakeFiles/psps_functions.dir/depend:
	cd /cluster/home/pafloria/testMoe5/buildO3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cluster/home/pafloria/testMoe5/Moe /cluster/home/pafloria/testMoe5/Moe /cluster/home/pafloria/testMoe5/buildO3 /cluster/home/pafloria/testMoe5/buildO3 /cluster/home/pafloria/testMoe5/buildO3/CMakeFiles/psps_functions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/psps_functions.dir/depend

