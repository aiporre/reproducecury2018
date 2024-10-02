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
CMAKE_SOURCE_DIR = /home/sauron/Documents/Phd/code/cury/catchine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp

# Include any dependencies generated for this target.
include CMakeFiles/testGPU_Cauchy3D.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/testGPU_Cauchy3D.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/testGPU_Cauchy3D.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testGPU_Cauchy3D.dir/flags.make

CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o: CMakeFiles/testGPU_Cauchy3D.dir/flags.make
CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o: ../src/testGPU_Cauchy.cpp
CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o: CMakeFiles/testGPU_Cauchy3D.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o -MF CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o.d -o CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o -c /home/sauron/Documents/Phd/code/cury/catchine/src/testGPU_Cauchy.cpp

CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sauron/Documents/Phd/code/cury/catchine/src/testGPU_Cauchy.cpp > CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.i

CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sauron/Documents/Phd/code/cury/catchine/src/testGPU_Cauchy.cpp -o CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.s

# Object files for target testGPU_Cauchy3D
testGPU_Cauchy3D_OBJECTS = \
"CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o"

# External object files for target testGPU_Cauchy3D
testGPU_Cauchy3D_EXTERNAL_OBJECTS =

../bin/testGPU_Cauchy3D: CMakeFiles/testGPU_Cauchy3D.dir/src/testGPU_Cauchy.cpp.o
../bin/testGPU_Cauchy3D: CMakeFiles/testGPU_Cauchy3D.dir/build.make
../bin/testGPU_Cauchy3D: /usr/local/cuda-11.8/lib64/libcudart_static.a
../bin/testGPU_Cauchy3D: /usr/lib/x86_64-linux-gnu/librt.a
../bin/testGPU_Cauchy3D: libCudaConvolution.a
../bin/testGPU_Cauchy3D: /usr/local/cuda-11.8/lib64/libcudart_static.a
../bin/testGPU_Cauchy3D: /usr/lib/x86_64-linux-gnu/librt.a
../bin/testGPU_Cauchy3D: CMakeFiles/testGPU_Cauchy3D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/testGPU_Cauchy3D"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testGPU_Cauchy3D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testGPU_Cauchy3D.dir/build: ../bin/testGPU_Cauchy3D
.PHONY : CMakeFiles/testGPU_Cauchy3D.dir/build

CMakeFiles/testGPU_Cauchy3D.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testGPU_Cauchy3D.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testGPU_Cauchy3D.dir/clean

CMakeFiles/testGPU_Cauchy3D.dir/depend:
	cd /home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sauron/Documents/Phd/code/cury/catchine /home/sauron/Documents/Phd/code/cury/catchine /home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp /home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp /home/sauron/Documents/Phd/code/cury/catchine/cmake_tmp/CMakeFiles/testGPU_Cauchy3D.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testGPU_Cauchy3D.dir/depend

