CC=clang++
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17
CFLAGS=-O3 --std=c++17
SOURCE_DIR=.
OBJECT_DIR=objects
DBG_NAME=debug
BIN_NAME=build
SIMPLE_PM_NAME=simple_PM


SIMPLE_PM_SOURCES=simple_PM LinearAlgebra Printers


SIMPLE_PM_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SIMPLE_PM_SOURCES)))
DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))

all: build
$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $@


simple_PM:
$(SIMPLE_PM_NAME): $(SIMPLE_PM_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out



debug:
$(DBG_NAME): $(DBG_OBJECTS)
	$(CC) $(CFLAGS_DBG) $^ -o $@.out

$(OBJECT_DIR)/%_DBG.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS_DBG) -c $< -o $@

clean:
	rm -f $(OBJECT_DIR)/*.o

purge:
	rm -rf $(OBJECT_DIR)/; rm $(SOURCE_DIR)/*.out
