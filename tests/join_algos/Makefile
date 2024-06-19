# Compiler
CXX := g++

# Directories
INCLUDE_DIR := ../../include
SRC_DIR := ../../src
JOIN_ALGOS_DIR := $(SRC_DIR)/join_algos
UTILS_DIR := $(SRC_DIR)/utils
TESTS_DIR := .
BUILD_DIR := build

# Flags
CXXFLAGS := -std=c++11 -I$(INCLUDE_DIR)

# Source and Object Files
JOIN_ALGOS_SRCS := $(wildcard $(JOIN_ALGOS_DIR)/*.cc)
UTILS_SRCS := $(wildcard $(UTILS_DIR)/*.cc)
TEST_SRCS := $(TESTS_DIR)/test_inter_filter.cc
SRCS := $(JOIN_ALGOS_SRCS) $(UTILS_SRCS)
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(BUILD_DIR)/%.o)
TEST_OBJS := $(TEST_SRCS:$(TESTS_DIR)/%.cc=$(BUILD_DIR)/%.o)

# Executable
TEST_EXEC := $(BUILD_DIR)/test_exec

# Rules
all: $(TEST_EXEC)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TESTS_DIR)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_EXEC): $(OBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

run: all
	./$(TEST_EXEC)

# Phony Targets
.PHONY: clean

clean:
	rm -rf $(BUILD_DIR) $(TEST_EXEC)