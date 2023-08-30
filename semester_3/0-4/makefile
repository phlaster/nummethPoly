CC=g++#clang++
CFLAGS_DBG=-g3 -O0 -Wall -Wextra --std=c++17
CFLAGS=-O3 --std=c++17
SOURCE_DIR=.
OBJECT_DIR=objects
DBG_NAME=debug
BIN_NAME=build
SAMPLE_NAME=sample
TASK_3_NAME=task3

SOURCES=main LinearAlgebra Printers
SAMPLE_SOURSES=report_sample LinearAlgebra Printers
TASK_3_SOURCES=task_3 LinearAlgebra Printers

TASK_3_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(TASK_3_SOURCES)))
SAMPLE_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SAMPLE_SOURSES)))
DBG_OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix _DBG.o,$(SOURCES)))
OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))


all: build
$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $@


sample:
$(SAMPLE_NAME): $(SAMPLE_OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out

task3:
$(TASK_3_NAME): $(TASK_3_OBJECTS)
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
