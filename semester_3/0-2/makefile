CC=clang++
CFLAGS=-g3
SOURCE_DIR=.
OBJECT_DIR=objects
BIN_NAME=build

SOURCES=main functions algorithms graphicBuild

OBJECTS=$(addprefix $(OBJECT_DIR)/,$(addsuffix .o,$(SOURCES)))

all: build

$(BIN_NAME): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@.out

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $@

clean:
	rm -f $(OBJECT_DIR)/*.o

purge:
	rm -rf $(OBJECT_DIR)/; rm $(SOURCE_DIR)/$(BIN_NAME).out