CC=g++
CXXFLAGS= -Wall -O2 -std=c++11
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LDFLAGS = -framework Accelerate
endif

ifeq ($(UNAME_S),Linux)
	LDFLAGS = -llapack -lblas
endif


all: main
	
main: main.o opt.o test_liag.o def.o dmatrix.o 
	$(CC) $(CFLAGS) $(LDFLAGS) main.o opt.o test_liag.o def.o dmatrix.o  -o main

main.o: main.cpp
	$(CC) -c main.cpp

opt.o: opt.cpp
	$(CC) -c opt.cpp
	
test_liag.o: test_liag.cpp
	$(CC) -c test_liag.cpp
	
def.o: def.c 
	$(CC) -c def.c

dmatrix.o: dmatrix.c
	$(CC) -c dmatrix.c

clean: 
	rm -vf *.o main
