# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /snap/clion/99/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/99/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/FaNDOM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FaNDOM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FaNDOM.dir/flags.make

CMakeFiles/FaNDOM.dir/main.cpp.o: CMakeFiles/FaNDOM.dir/flags.make
CMakeFiles/FaNDOM.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/FaNDOM.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FaNDOM.dir/main.cpp.o -c /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/main.cpp

CMakeFiles/FaNDOM.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FaNDOM.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/main.cpp > CMakeFiles/FaNDOM.dir/main.cpp.i

CMakeFiles/FaNDOM.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FaNDOM.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/main.cpp -o CMakeFiles/FaNDOM.dir/main.cpp.s

CMakeFiles/FaNDOM.dir/OMIO.cpp.o: CMakeFiles/FaNDOM.dir/flags.make
CMakeFiles/FaNDOM.dir/OMIO.cpp.o: ../OMIO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/FaNDOM.dir/OMIO.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FaNDOM.dir/OMIO.cpp.o -c /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/OMIO.cpp

CMakeFiles/FaNDOM.dir/OMIO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FaNDOM.dir/OMIO.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/OMIO.cpp > CMakeFiles/FaNDOM.dir/OMIO.cpp.i

CMakeFiles/FaNDOM.dir/OMIO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FaNDOM.dir/OMIO.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/OMIO.cpp -o CMakeFiles/FaNDOM.dir/OMIO.cpp.s

CMakeFiles/FaNDOM.dir/OMHelper.cpp.o: CMakeFiles/FaNDOM.dir/flags.make
CMakeFiles/FaNDOM.dir/OMHelper.cpp.o: ../OMHelper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/FaNDOM.dir/OMHelper.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FaNDOM.dir/OMHelper.cpp.o -c /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/OMHelper.cpp

CMakeFiles/FaNDOM.dir/OMHelper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FaNDOM.dir/OMHelper.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/OMHelper.cpp > CMakeFiles/FaNDOM.dir/OMHelper.cpp.i

CMakeFiles/FaNDOM.dir/OMHelper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FaNDOM.dir/OMHelper.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/OMHelper.cpp -o CMakeFiles/FaNDOM.dir/OMHelper.cpp.s

CMakeFiles/FaNDOM.dir/fandomAligner.cpp.o: CMakeFiles/FaNDOM.dir/flags.make
CMakeFiles/FaNDOM.dir/fandomAligner.cpp.o: ../fandomAligner.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/FaNDOM.dir/fandomAligner.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FaNDOM.dir/fandomAligner.cpp.o -c /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/fandomAligner.cpp

CMakeFiles/FaNDOM.dir/fandomAligner.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FaNDOM.dir/fandomAligner.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/fandomAligner.cpp > CMakeFiles/FaNDOM.dir/fandomAligner.cpp.i

CMakeFiles/FaNDOM.dir/fandomAligner.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FaNDOM.dir/fandomAligner.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/fandomAligner.cpp -o CMakeFiles/FaNDOM.dir/fandomAligner.cpp.s

# Object files for target FaNDOM
FaNDOM_OBJECTS = \
"CMakeFiles/FaNDOM.dir/main.cpp.o" \
"CMakeFiles/FaNDOM.dir/OMIO.cpp.o" \
"CMakeFiles/FaNDOM.dir/OMHelper.cpp.o" \
"CMakeFiles/FaNDOM.dir/fandomAligner.cpp.o"

# External object files for target FaNDOM
FaNDOM_EXTERNAL_OBJECTS =

FaNDOM: CMakeFiles/FaNDOM.dir/main.cpp.o
FaNDOM: CMakeFiles/FaNDOM.dir/OMIO.cpp.o
FaNDOM: CMakeFiles/FaNDOM.dir/OMHelper.cpp.o
FaNDOM: CMakeFiles/FaNDOM.dir/fandomAligner.cpp.o
FaNDOM: CMakeFiles/FaNDOM.dir/build.make
FaNDOM: CMakeFiles/FaNDOM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable FaNDOM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FaNDOM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FaNDOM.dir/build: FaNDOM

.PHONY : CMakeFiles/FaNDOM.dir/build

CMakeFiles/FaNDOM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FaNDOM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FaNDOM.dir/clean

CMakeFiles/FaNDOM.dir/depend:
	cd /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug /home/jens/Dropbox/BafnaLab/bionano_analysis/method_development/FaNDOM/cmake-build-debug/CMakeFiles/FaNDOM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FaNDOM.dir/depend

