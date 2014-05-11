TARGET  = bst_driver
USERCFG = config.make
SRC     = $(wildcard opt/bst_*.c) driver.c driver_log.c ref/bst_ref.c urange.c 
OBJ     = $(patsubst %.c, %.o, $(SRC))

include $(USERCFG) # for DEFINES and FLAGS

CC        = gcc
LIBS      = -lm -lpfm
INC       = -I.
USERFLAGS = $(DEFINES) $(FLAGS)
#GITREV    = $(shell git rev-parse HEAD)
GITREV    = $(shell git log --oneline -n1)
M_ENV     = -DM_ENV_USERFLAGS="$(USERFLAGS)" -DM_ENV_GITREV="$(GITREV)"
CF        = -g -std=c99 $(INC) $(USERFLAGS) $(M_ENV)

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
