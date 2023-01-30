CC = gcc

IPATH = -I/usr/local/include -L/usr/local/lib

LIBS = -lm -lgsl -lgslcblas

CFLAGS = -g -Wall -O2 -funroll-loops $(IPATH)

all: cdt

cdt: CDT2d.c main.c
	$(CC) $(CFLAGS) CDT2d.c main.c $(LIBS) -o cdt
