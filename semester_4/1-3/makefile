CC=clang++
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17
CFLAGS=-O3 --std=c++17
SOURCE_DIR=.
OBJECT_DIR=objects
CSV_DIR=CSVs
DBG_NAME=debug
BIN_NAME=build
SAMPLE_NAME=sample
MULTI_NAME=multi
TESTS_NAME=tests


SOURCES=main Functions Printers PCG
SAMPLE_SOURCES=sample Functions Printers PCG
MULTI_SOURCES=multi Functions Printers PCG
TESTS_SOURCES=tests Functions Printers PCG

DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))
SAMPLE_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SAMPLE_SOURCES)))
MULTI_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(MULTI_SOURCES)))
TESTS_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(TESTS_SOURCES)))


all: build

$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; mkdir $(CSV_DIR)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	if ! test -d $(OBJECT_DIR); then mkdir -p $(OBJECT_DIR); fi


sample:
$(SAMPLE_NAME): $(SAMPLE_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; if ! test -d $(CSV_DIR); then mkdir -p $(CSV_DIR); fi
$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@


multi:
$(MULTI_NAME): $(MULTI_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; if ! test -d $(CSV_DIR); then mkdir -p $(CSV_DIR); fi
$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

tests:
$(TESTS_NAME): $(TESTS_OBJECTS)
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