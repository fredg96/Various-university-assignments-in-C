###############################################################################
# Makefile for individual project, Parallel and Distributed Computing 2020.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = malaria

all: $(BIN)

malaria: malaria.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
