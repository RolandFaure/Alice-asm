# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_SOURCE_DIR = /home/rfaure/Documents/these/MSR/wonderasm/bcalm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build

# Utility rule file for thirdparty_copy.

# Include any custom commands dependencies for this target.
include ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/compiler_depend.make

# Include the progress variables for this target.
include ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/progress.make

thirdparty_copy: ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/build.make
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E make_directory /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//json
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E make_directory /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//BooPHF
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E make_directory /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//kff-cpp-api
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E copy_if_different /home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/BooPHF/BooPHF.h /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//BooPHF/
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E copy_if_different /home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/json/json.hpp /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//json/
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E copy_if_different /home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/kff-cpp-api/kff_io.hpp /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//kff-cpp-api/
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && /usr/bin/cmake -E copy_if_different /home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/kff-cpp-api/merge.hpp /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include//kff-cpp-api/
.PHONY : thirdparty_copy

# Rule to build all files generated by this target.
ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/build: thirdparty_copy
.PHONY : ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/build

ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/clean:
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty && $(CMAKE_COMMAND) -P CMakeFiles/thirdparty_copy.dir/cmake_clean.cmake
.PHONY : ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/clean

ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/depend:
	cd /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rfaure/Documents/these/MSR/wonderasm/bcalm /home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty /home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/depend

