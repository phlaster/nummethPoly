CC=clang++
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17
CFLAGS=-O3 --std=c++17
SOURCE_DIR=.
OBJECT_DIR=objects
CSV_DIR=CSVs
DBG_NAME=debug
BIN_NAME=build

SOURCES=main functions algorithms graphicBuild printers Buffer tasks

DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))

all: build

$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; mkdir $(CSV_DIR)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $@



debug:

$(DBG_NAME): $(DBG_OBJECTS)
	$(CC) $(CFLAGS_DBG) $^ -o $@.out; mkdir $(CSV_DIR)

$(OBJECT_DIR)/%_DBG.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS_DBG) -c $< -o $@




clean:
	rm -f $(OBJECT_DIR)/*.o

purge:
	rm -rf $(OBJECT_DIR)/; rm -rf $(CSV_DIR); rm $(SOURCE_DIR)/*.out