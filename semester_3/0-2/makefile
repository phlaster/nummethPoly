CC=g++ #clang++
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17
CFLAGS=-O3 --std=c++17
SOURCE_DIR=.
OBJECT_DIR=objects
CSV_DIR=CSVs
DBG_NAME=debug
BIN_NAME=build
SINGLEDOT_NAME=singledot

SOURCES=main functions algorithms graphicBuild printers Buffer tasks
SOURCES_SINGLEDOT=track_single_dot functions algorithms graphicBuild printers Buffer tasks

DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))
SINGLEDOT_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES_SINGLEDOT)))

all: build

$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; mkdir $(CSV_DIR)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	if ! test -d $(OBJECT_DIR); then mkdir -p $(OBJECT_DIR); fi


singledot:
$(SINGLEDOT_NAME): $(SINGLEDOT_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; if ! test -d $(CSV_DIR); then mkdir -p $(CSV_DIR); fi
$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@



debug:
$(DBG_NAME): $(DBG_OBJECTS)
	$(CC) $(CFLAGS_DBG) $^ -o $@.out; if ! test -d $(CSV_DIR); then mkdir -p $(CSV_DIR); fi
$(OBJECT_DIR)/%_DBG.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS_DBG) -c $< -o $@




clean:
	rm -f $(OBJECT_DIR)/*.o

purge:
	rm -rf $(OBJECT_DIR)/; rm -rf $(CSV_DIR); rm $(SOURCE_DIR)/*.out