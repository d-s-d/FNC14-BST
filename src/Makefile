TARGET  = bst_driver
USERCFG = config.make
SRC     = driver.c ref/bst_ref.c
OBJ     = $(patsubst %.c, %.o, $(SRC))

include $(USERCFG) # for DEFINES and FLAGS

CC        = gcc
LIBS      = -lm -lpfm
INC       = -I.
USERFLAGS = $(DEFINES) $(FLAGS)
CF        = -g -std=c99 $(INC) $(USERFLAGS)

.PHONY: all clean

all: $(TARGET)

clean:
	rm -f $(OBJ)
	rm -f $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

%.o: %.c $(USERCFG)
	$(CC) $(CF) -c $< -o $@

$(USERCFG):
	cp $(USERCFG).example $(USERCFG)
