# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build

# Include any dependencies generated for this target.
include ConvexityTests/CMakeFiles/ConvexityTests.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include ConvexityTests/CMakeFiles/ConvexityTests.dir/compiler_depend.make

# Include the progress variables for this target.
include ConvexityTests/CMakeFiles/ConvexityTests.dir/progress.make

# Include the compile flags for this target's objects.
include ConvexityTests/CMakeFiles/ConvexityTests.dir/flags.make

ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.o: ConvexityTests/CMakeFiles/ConvexityTests.dir/flags.make
ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.o: /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/ConvexityTests/main.cc
ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.o: ConvexityTests/CMakeFiles/ConvexityTests.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.o"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.o -MF CMakeFiles/ConvexityTests.dir/main.cc.o.d -o CMakeFiles/ConvexityTests.dir/main.cc.o -c /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/ConvexityTests/main.cc

ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/ConvexityTests.dir/main.cc.i"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/ConvexityTests/main.cc > CMakeFiles/ConvexityTests.dir/main.cc.i

ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/ConvexityTests.dir/main.cc.s"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/ConvexityTests/main.cc -o CMakeFiles/ConvexityTests.dir/main.cc.s

# Object files for target ConvexityTests
ConvexityTests_OBJECTS = \
"CMakeFiles/ConvexityTests.dir/main.cc.o"

# External object files for target ConvexityTests
ConvexityTests_EXTERNAL_OBJECTS =

Build/bin/ConvexityTests: ConvexityTests/CMakeFiles/ConvexityTests.dir/main.cc.o
Build/bin/ConvexityTests: ConvexityTests/CMakeFiles/ConvexityTests.dir/build.make
Build/bin/ConvexityTests: ConvexityTests/CMakeFiles/ConvexityTests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/ConvexityTests"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ConvexityTests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ConvexityTests/CMakeFiles/ConvexityTests.dir/build: Build/bin/ConvexityTests
.PHONY : ConvexityTests/CMakeFiles/ConvexityTests.dir/build

ConvexityTests/CMakeFiles/ConvexityTests.dir/clean:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests && $(CMAKE_COMMAND) -P CMakeFiles/ConvexityTests.dir/cmake_clean.cmake
.PHONY : ConvexityTests/CMakeFiles/ConvexityTests.dir/clean

ConvexityTests/CMakeFiles/ConvexityTests.dir/depend:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/ConvexityTests /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/ConvexityTests/CMakeFiles/ConvexityTests.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : ConvexityTests/CMakeFiles/ConvexityTests.dir/depend

