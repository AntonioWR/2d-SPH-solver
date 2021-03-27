lib=-lboost_program_options
head=SPH.h

default: myprog

main.o:main.cpp
	mpicxx -Wall -O2 -g -o main.o -c main.cpp

myclass.o:SPH.cpp $(head)
	mpicxx -Wall -O2 -g -o myclass.o -c SPH.cpp

myprog:main.o myclass.o $(head)
	mpicxx -o myprog main.o myclass.o $(lib)


.PHONY: clean run


