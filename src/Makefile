TARGET = bst_driver
SRC    = driver.c ref/bst_ref.c
OBJ    = $(patsubst %.c, %.o, $(SRC))

CC     = gcc
LIBS   = -lm -lpfm
INC    = -I.
CF     = -g -std=c99 $(INC)

.PHONY: all clean

all: $(TARGET)

clean:
	rm -f $(OBJ)
	rm -f $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

%.o: %.c
	$(CC) $(CF) -c $< -o $@
