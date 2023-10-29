CC=clang++

# Flags for build
CFLAGS=-O3 --std=c++17

# Flags for debugging
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17

# Where to take sources
SOURCE_DIR=.

# Where to put objects
OBJECT_DIR=objects

DBG_NAME=debug
BIN_NAME=build

# Names of your bins
SOURCES=main functions routines runge

# Baking names for objects for build & objects for debug
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))
DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))


# Everything below is kinda invariant, you can
all: build
$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out
$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJECT_DIR):
	mkdir -p $@

debug:
$(DBG_NAME): $(DBG_OBJECTS)
	$(CC) $(CFLAGS_DBG) $^ -o $@.out
$(OBJECT_DIR)/%_DBG.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS_DBG) -c $< -o $@

clean: # only deleting the object files
	rm -f $(OBJECT_DIR)/*.o

purge:  # all csvs will be deleted in cwd!
	rm -rf $(OBJECT_DIR)/; rm $(SOURCE_DIR)/*.out; rm $(SOURCE_DIR)/*.csv
