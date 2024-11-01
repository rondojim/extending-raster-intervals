# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_SOURCE_DIR = /home/dimitriskpl/DSIT/2_nd_semester/DB_Systems/Project/phase1-paper-presentation-roadmap-2024-raster-intervals

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dimitriskpl/DSIT/2_nd_semester/DB_Systems/Project/phase1-paper-presentation-roadmap-2024-raster-intervals

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test
.PHONY : test/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/dimitriskpl/DSIT/2_nd_semester/DB_Systems/Project/phase1-paper-presentation-roadmap-2024-raster-intervals/CMakeFiles /home/dimitriskpl/DSIT/2_nd_semester/DB_Systems/Project/phase1-paper-presentation-roadmap-2024-raster-intervals//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/dimitriskpl/DSIT/2_nd_semester/DB_Systems/Project/phase1-paper-presentation-roadmap-2024-raster-intervals/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named main

# Build rule for target.
main: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 main
.PHONY : main

# fast build rule for target.
main/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/build
.PHONY : main/fast

src/join_algos/combined_sweep.o: src/join_algos/combined_sweep.cc.o
.PHONY : src/join_algos/combined_sweep.o

# target to build an object file
src/join_algos/combined_sweep.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/combined_sweep.cc.o
.PHONY : src/join_algos/combined_sweep.cc.o

src/join_algos/combined_sweep.i: src/join_algos/combined_sweep.cc.i
.PHONY : src/join_algos/combined_sweep.i

# target to preprocess a source file
src/join_algos/combined_sweep.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/combined_sweep.cc.i
.PHONY : src/join_algos/combined_sweep.cc.i

src/join_algos/combined_sweep.s: src/join_algos/combined_sweep.cc.s
.PHONY : src/join_algos/combined_sweep.s

# target to generate assembly for a file
src/join_algos/combined_sweep.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/combined_sweep.cc.s
.PHONY : src/join_algos/combined_sweep.cc.s

src/join_algos/inter_filter.o: src/join_algos/inter_filter.cc.o
.PHONY : src/join_algos/inter_filter.o

# target to build an object file
src/join_algos/inter_filter.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/inter_filter.cc.o
.PHONY : src/join_algos/inter_filter.cc.o

src/join_algos/inter_filter.i: src/join_algos/inter_filter.cc.i
.PHONY : src/join_algos/inter_filter.i

# target to preprocess a source file
src/join_algos/inter_filter.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/inter_filter.cc.i
.PHONY : src/join_algos/inter_filter.cc.i

src/join_algos/inter_filter.s: src/join_algos/inter_filter.cc.s
.PHONY : src/join_algos/inter_filter.s

# target to generate assembly for a file
src/join_algos/inter_filter.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/inter_filter.cc.s
.PHONY : src/join_algos/inter_filter.cc.s

src/join_algos/mbr_no_inside.o: src/join_algos/mbr_no_inside.cc.o
.PHONY : src/join_algos/mbr_no_inside.o

# target to build an object file
src/join_algos/mbr_no_inside.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_no_inside.cc.o
.PHONY : src/join_algos/mbr_no_inside.cc.o

src/join_algos/mbr_no_inside.i: src/join_algos/mbr_no_inside.cc.i
.PHONY : src/join_algos/mbr_no_inside.i

# target to preprocess a source file
src/join_algos/mbr_no_inside.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_no_inside.cc.i
.PHONY : src/join_algos/mbr_no_inside.cc.i

src/join_algos/mbr_no_inside.s: src/join_algos/mbr_no_inside.cc.s
.PHONY : src/join_algos/mbr_no_inside.s

# target to generate assembly for a file
src/join_algos/mbr_no_inside.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_no_inside.cc.s
.PHONY : src/join_algos/mbr_no_inside.cc.s

src/join_algos/mbr_no_outside.o: src/join_algos/mbr_no_outside.cc.o
.PHONY : src/join_algos/mbr_no_outside.o

# target to build an object file
src/join_algos/mbr_no_outside.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_no_outside.cc.o
.PHONY : src/join_algos/mbr_no_outside.cc.o

src/join_algos/mbr_no_outside.i: src/join_algos/mbr_no_outside.cc.i
.PHONY : src/join_algos/mbr_no_outside.i

# target to preprocess a source file
src/join_algos/mbr_no_outside.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_no_outside.cc.i
.PHONY : src/join_algos/mbr_no_outside.cc.i

src/join_algos/mbr_no_outside.s: src/join_algos/mbr_no_outside.cc.s
.PHONY : src/join_algos/mbr_no_outside.s

# target to generate assembly for a file
src/join_algos/mbr_no_outside.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_no_outside.cc.s
.PHONY : src/join_algos/mbr_no_outside.cc.s

src/join_algos/mbr_suboptimal_2d_segt.o: src/join_algos/mbr_suboptimal_2d_segt.cc.o
.PHONY : src/join_algos/mbr_suboptimal_2d_segt.o

# target to build an object file
src/join_algos/mbr_suboptimal_2d_segt.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_suboptimal_2d_segt.cc.o
.PHONY : src/join_algos/mbr_suboptimal_2d_segt.cc.o

src/join_algos/mbr_suboptimal_2d_segt.i: src/join_algos/mbr_suboptimal_2d_segt.cc.i
.PHONY : src/join_algos/mbr_suboptimal_2d_segt.i

# target to preprocess a source file
src/join_algos/mbr_suboptimal_2d_segt.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_suboptimal_2d_segt.cc.i
.PHONY : src/join_algos/mbr_suboptimal_2d_segt.cc.i

src/join_algos/mbr_suboptimal_2d_segt.s: src/join_algos/mbr_suboptimal_2d_segt.cc.s
.PHONY : src/join_algos/mbr_suboptimal_2d_segt.s

# target to generate assembly for a file
src/join_algos/mbr_suboptimal_2d_segt.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_suboptimal_2d_segt.cc.s
.PHONY : src/join_algos/mbr_suboptimal_2d_segt.cc.s

src/join_algos/mbr_sweep_brute.o: src/join_algos/mbr_sweep_brute.cc.o
.PHONY : src/join_algos/mbr_sweep_brute.o

# target to build an object file
src/join_algos/mbr_sweep_brute.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_sweep_brute.cc.o
.PHONY : src/join_algos/mbr_sweep_brute.cc.o

src/join_algos/mbr_sweep_brute.i: src/join_algos/mbr_sweep_brute.cc.i
.PHONY : src/join_algos/mbr_sweep_brute.i

# target to preprocess a source file
src/join_algos/mbr_sweep_brute.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_sweep_brute.cc.i
.PHONY : src/join_algos/mbr_sweep_brute.cc.i

src/join_algos/mbr_sweep_brute.s: src/join_algos/mbr_sweep_brute.cc.s
.PHONY : src/join_algos/mbr_sweep_brute.s

# target to generate assembly for a file
src/join_algos/mbr_sweep_brute.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/mbr_sweep_brute.cc.s
.PHONY : src/join_algos/mbr_sweep_brute.cc.s

src/join_algos/node.o: src/join_algos/node.cc.o
.PHONY : src/join_algos/node.o

# target to build an object file
src/join_algos/node.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/node.cc.o
.PHONY : src/join_algos/node.cc.o

src/join_algos/node.i: src/join_algos/node.cc.i
.PHONY : src/join_algos/node.i

# target to preprocess a source file
src/join_algos/node.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/node.cc.i
.PHONY : src/join_algos/node.cc.i

src/join_algos/node.s: src/join_algos/node.cc.s
.PHONY : src/join_algos/node.s

# target to generate assembly for a file
src/join_algos/node.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/join_algos/node.cc.s
.PHONY : src/join_algos/node.cc.s

src/main.o: src/main.cc.o
.PHONY : src/main.o

# target to build an object file
src/main.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/main.cc.o
.PHONY : src/main.cc.o

src/main.i: src/main.cc.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/main.cc.i
.PHONY : src/main.cc.i

src/main.s: src/main.cc.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/main.cc.s
.PHONY : src/main.cc.s

src/utils/data_reader.o: src/utils/data_reader.cc.o
.PHONY : src/utils/data_reader.o

# target to build an object file
src/utils/data_reader.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/utils/data_reader.cc.o
.PHONY : src/utils/data_reader.cc.o

src/utils/data_reader.i: src/utils/data_reader.cc.i
.PHONY : src/utils/data_reader.i

# target to preprocess a source file
src/utils/data_reader.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/utils/data_reader.cc.i
.PHONY : src/utils/data_reader.cc.i

src/utils/data_reader.s: src/utils/data_reader.cc.s
.PHONY : src/utils/data_reader.s

# target to generate assembly for a file
src/utils/data_reader.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/utils/data_reader.cc.s
.PHONY : src/utils/data_reader.cc.s

src/utils/geometry_types.o: src/utils/geometry_types.cc.o
.PHONY : src/utils/geometry_types.o

# target to build an object file
src/utils/geometry_types.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/utils/geometry_types.cc.o
.PHONY : src/utils/geometry_types.cc.o

src/utils/geometry_types.i: src/utils/geometry_types.cc.i
.PHONY : src/utils/geometry_types.i

# target to preprocess a source file
src/utils/geometry_types.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/utils/geometry_types.cc.i
.PHONY : src/utils/geometry_types.cc.i

src/utils/geometry_types.s: src/utils/geometry_types.cc.s
.PHONY : src/utils/geometry_types.s

# target to generate assembly for a file
src/utils/geometry_types.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/utils/geometry_types.cc.s
.PHONY : src/utils/geometry_types.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... test"
	@echo "... main"
	@echo "... src/join_algos/combined_sweep.o"
	@echo "... src/join_algos/combined_sweep.i"
	@echo "... src/join_algos/combined_sweep.s"
	@echo "... src/join_algos/inter_filter.o"
	@echo "... src/join_algos/inter_filter.i"
	@echo "... src/join_algos/inter_filter.s"
	@echo "... src/join_algos/mbr_no_inside.o"
	@echo "... src/join_algos/mbr_no_inside.i"
	@echo "... src/join_algos/mbr_no_inside.s"
	@echo "... src/join_algos/mbr_no_outside.o"
	@echo "... src/join_algos/mbr_no_outside.i"
	@echo "... src/join_algos/mbr_no_outside.s"
	@echo "... src/join_algos/mbr_suboptimal_2d_segt.o"
	@echo "... src/join_algos/mbr_suboptimal_2d_segt.i"
	@echo "... src/join_algos/mbr_suboptimal_2d_segt.s"
	@echo "... src/join_algos/mbr_sweep_brute.o"
	@echo "... src/join_algos/mbr_sweep_brute.i"
	@echo "... src/join_algos/mbr_sweep_brute.s"
	@echo "... src/join_algos/node.o"
	@echo "... src/join_algos/node.i"
	@echo "... src/join_algos/node.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/utils/data_reader.o"
	@echo "... src/utils/data_reader.i"
	@echo "... src/utils/data_reader.s"
	@echo "... src/utils/geometry_types.o"
	@echo "... src/utils/geometry_types.i"
	@echo "... src/utils/geometry_types.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

