#CC=gcc -DM_FLOPS -DM_MISSES
#CC=gcc -DM_FLOPS
#CC=gcc -DM_MISSES
CC=gcc -g
MARCH=native

all: mmm.o main.o
	$(CC) -o mmmperf mmm.o main.o

mmm.o: mmm.c
	$(CC) -march=$(MARCH) -c -std=c99 -O3 -o mmm.o mmm.c

main.o: main.c
	$(CC) -c -o main.o main.c
	
clean:
	rm -f *.o mmmperf
