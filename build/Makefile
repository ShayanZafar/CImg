# Makefile for w1
#
GCC_VERSION = 4.7.1
PREFIX = /usr/local/gcc/${GCC_VERSION}/bin/
CC = ${PREFIX}gcc
CPP = ${PREFIX}g++

triangle: triangle.o
	$(CPP) -pg -ow1 triangle.o

triangle.o: triangle.cpp
	$(CPP) -c -O2 -g -pg -std=c++0x triangle.cpp

clean:
	rm *.o
