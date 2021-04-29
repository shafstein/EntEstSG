# Makefile adapted from https://spin.atomicobject.com/2016/08/26/makefile-c-projects/ (Job Vranish)
TARGET_EXEC ?= EstEntSGex
CC = g++
BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
CXXFLAGS = -O2 
LFLAGS =-larmadillo -pthread
SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

$(TARGET_EXEC): $(OBJS)
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LFLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ $(LFLAGS)


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p

