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
CMAKE_BINARY_DIR = /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11

# Include any dependencies generated for this target.
include EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/compiler_depend.make

# Include the progress variables for this target.
include EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/progress.make

# Include the compile flags for this target's objects.
include EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/flags.make

EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o: EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/flags.make
EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o: EqualityConstrainedNewton/unit_tests.cc
EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o: EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o -MF CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o.d -o CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o -c /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton/unit_tests.cc

EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.i"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton/unit_tests.cc > CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.i

EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.s"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton/unit_tests.cc -o CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.s

# Object files for target EqualityConstrainedNewton-test
EqualityConstrainedNewton__test_OBJECTS = \
"CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o"

# External object files for target EqualityConstrainedNewton-test
EqualityConstrainedNewton__test_EXTERNAL_OBJECTS =

Build/bin/EqualityConstrainedNewton-test: EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/unit_tests.cc.o
Build/bin/EqualityConstrainedNewton-test: EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/build.make
Build/bin/EqualityConstrainedNewton-test: lib/libgtest.a
Build/bin/EqualityConstrainedNewton-test: lib/libgtest_main.a
Build/bin/EqualityConstrainedNewton-test: lib/libgtest.a
Build/bin/EqualityConstrainedNewton-test: EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/EqualityConstrainedNewton-test"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/EqualityConstrainedNewton-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/build: Build/bin/EqualityConstrainedNewton-test
.PHONY : EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/build

EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/clean:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton && $(CMAKE_COMMAND) -P CMakeFiles/EqualityConstrainedNewton-test.dir/cmake_clean.cmake
.PHONY : EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/clean

EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/depend:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : EqualityConstrainedNewton/CMakeFiles/EqualityConstrainedNewton-test.dir/depend

