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
include NewtonMethods/CMakeFiles/NewtonMethods.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include NewtonMethods/CMakeFiles/NewtonMethods.dir/compiler_depend.make

# Include the progress variables for this target.
include NewtonMethods/CMakeFiles/NewtonMethods.dir/progress.make

# Include the compile flags for this target's objects.
include NewtonMethods/CMakeFiles/NewtonMethods.dir/flags.make

NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.o: NewtonMethods/CMakeFiles/NewtonMethods.dir/flags.make
NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.o: NewtonMethods/main.cc
NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.o: NewtonMethods/CMakeFiles/NewtonMethods.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.o"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.o -MF CMakeFiles/NewtonMethods.dir/main.cc.o.d -o CMakeFiles/NewtonMethods.dir/main.cc.o -c /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/main.cc

NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/NewtonMethods.dir/main.cc.i"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/main.cc > CMakeFiles/NewtonMethods.dir/main.cc.i

NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/NewtonMethods.dir/main.cc.s"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/main.cc -o CMakeFiles/NewtonMethods.dir/main.cc.s

# Object files for target NewtonMethods
NewtonMethods_OBJECTS = \
"CMakeFiles/NewtonMethods.dir/main.cc.o"

# External object files for target NewtonMethods
NewtonMethods_EXTERNAL_OBJECTS =

Build/bin/NewtonMethods: NewtonMethods/CMakeFiles/NewtonMethods.dir/main.cc.o
Build/bin/NewtonMethods: NewtonMethods/CMakeFiles/NewtonMethods.dir/build.make
Build/bin/NewtonMethods: NewtonMethods/CMakeFiles/NewtonMethods.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../Build/bin/NewtonMethods"
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NewtonMethods.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
NewtonMethods/CMakeFiles/NewtonMethods.dir/build: Build/bin/NewtonMethods
.PHONY : NewtonMethods/CMakeFiles/NewtonMethods.dir/build

NewtonMethods/CMakeFiles/NewtonMethods.dir/clean:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods && $(CMAKE_COMMAND) -P CMakeFiles/NewtonMethods.dir/cmake_clean.cmake
.PHONY : NewtonMethods/CMakeFiles/NewtonMethods.dir/clean

NewtonMethods/CMakeFiles/NewtonMethods.dir/depend:
	cd /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11 /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeFiles/NewtonMethods.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : NewtonMethods/CMakeFiles/NewtonMethods.dir/depend

