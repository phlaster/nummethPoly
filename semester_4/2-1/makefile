CC=clang++
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17
CFLAGS=-O3 --std=c++17
SOURCE_DIR=.
OBJECT_DIR=objects
DBG_NAME=debug
MULTI_NAME=multi
BIN_NAME=build


SOURCES=main LinearAlgebra Printers algos
SOURCES_MULTI=multi LinearAlgebra Printers algos

MULTI_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES_MULTI)))
DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))

all: build
$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $@

multi:
$(MULTI_NAME): $(MULTI_OBJECTS)
	$(CC) $(CFLAGS_DBG) $^ -o $@.out

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS_DBG) -c $< -o $@


debug:
$(DBG_NAME): $(DBG_OBJECTS)
	$(CC) $(CFLAGS_DBG) $^ -o $@.out

$(OBJECT_DIR)/%_DBG.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS_DBG) -c $< -o $@

clean:
	rm -f $(OBJECT_DIR)/*.o

purge:
	rm -rf $(OBJECT_DIR)/; rm $(SOURCE_DIR)/*.out
