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
include GradientDescent/CMakeFiles/GradientDescent-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include GradientDescent/CMakeFiles/GradientDescent-test.dir/compiler_depend.make

# Include the progress variables for this target.
include GradientDescent/CMakeFiles/GradientDescent-test.dir/progress.make

# Include the compile flags for this target's objects.
include GradientDescent/CMakeFiles/GradientDescent-test.dir/flags.make

GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o: GradientDescent/CMakeFiles/GradientDescent-test.dir/flags.make
GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o: GradientDescent/unit_tests.cc
GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o: GradientDescent/CMakeFiles/GradientDescent-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o -MF CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o.d -o CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o -c /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent/unit_tests.cc

GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/GradientDescent-test.dir/unit_tests.cc.i"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent/unit_tests.cc > CMakeFiles/GradientDescent-test.dir/unit_tests.cc.i

GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/GradientDescent-test.dir/unit_tests.cc.s"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent/unit_tests.cc -o CMakeFiles/GradientDescent-test.dir/unit_tests.cc.s

# Object files for target GradientDescent-test
GradientDescent__test_OBJECTS = \
"CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o"

# External object files for target GradientDescent-test
GradientDescent__test_EXTERNAL_OBJECTS =

Build/bin/GradientDescent-test: GradientDescent/CMakeFiles/GradientDescent-test.dir/unit_tests.cc.o
Build/bin/GradientDescent-test: GradientDescent/CMakeFiles/GradientDescent-test.dir/build.make
Build/bin/GradientDescent-test: lib/libgtest.a
Build/bin/GradientDescent-test: lib/libgtest_main.a
Build/bin/GradientDescent-test: lib/libgtest.a
Build/bin/GradientDescent-test: GradientDescent/CMakeFiles/GradientDescent-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/GradientDescent-test"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GradientDescent-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
GradientDescent/CMakeFiles/GradientDescent-test.dir/build: Build/bin/GradientDescent-test
.PHONY : GradientDescent/CMakeFiles/GradientDescent-test.dir/build

GradientDescent/CMakeFiles/GradientDescent-test.dir/clean:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent && $(CMAKE_COMMAND) -P CMakeFiles/GradientDescent-test.dir/cmake_clean.cmake
.PHONY : GradientDescent/CMakeFiles/GradientDescent-test.dir/clean

GradientDescent/CMakeFiles/GradientDescent-test.dir/depend:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/GradientDescent/CMakeFiles/GradientDescent-test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : GradientDescent/CMakeFiles/GradientDescent-test.dir/depend

