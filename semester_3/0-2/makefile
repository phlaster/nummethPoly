CC=clang++
CFLAGS=-g3
SOURCE_DIR=.
OBJECT_DIR=objects
CSV_DIR=CSVs
BIN_NAME=build

SOURCES=main functions algorithms graphicBuild printers Buffer

OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))

all: build

$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out; mkdir $(CSV_DIR)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $@

clean:
	rm -f $(OBJECT_DIR)/*.o

purge:
	rm -rf $(OBJECT_DIR)/; rm -rf $(CSV_DIR); rm $(SOURCE_DIR)/$(BIN_NAME).out