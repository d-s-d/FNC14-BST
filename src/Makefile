TARGET  = bst_driver
RL_TARGET = roofline_driver
USERCFG = config.make
SRC     = $(wildcard opt/bst_*.c) driver.c driver_log.c ref/bst_ref.c urange.c
SRC_RL  = $(wildcard opt/bst_*.c) roofline_driver.c driver_log.c ref/bst_ref.c urange.c
OBJ     = $(patsubst %.c, %.o, $(SRC))
OBJ_RL  = $(patsubst %.c, %.o, $(SRC_RL))

include $(USERCFG) # for DEFINES and FLAGS

CC        = gcc
LIBS      = -lm -lpfm -pthread
INC       = -I.
USERFLAGS = $(DEFINES) $(FLAGS)
#GITREV    = $(shell git rev-parse HEAD)
GITREV    = $(shell git log --oneline -n1)
M_ENV     = -DM_ENV_USERFLAGS="\"$(USERFLAGS)\"" -DM_ENV_GITREV="\"$(GITREV)\""
CF        = -g -std=c99 $(INC) $(USERFLAGS) $(M_ENV)

.PHONY: all clean

all: $(TARGET) $(RL_TARGET)

clean:
	rm -f $(OBJ)
	rm -f $(OBJ_RL)
	rm -f $(TARGET)
	rm -f $(RL_TARGET)

$(TARGET): $(OBJ) variants.h
	$(CC) $^  $(LIBS) -o $@

$(RL_TARGET): $(OBJ_RL) variants.h MeasuringCore.lib
	g++ $^ $(LIBS) -o $@

%.o: %.c $(USERCFG)
	$(CC) $(CF) -c $< -o $@

$(USERCFG):
	cp $(USERCFG).example $(USERCFG)
