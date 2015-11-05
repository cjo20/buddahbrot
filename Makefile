OBJS = main.cpp
OBJ_NAME = main
CC = g++
COMPILE_FLAGS = -w -fopenmp -std=c++11 -O2
LINKER_FLAGS = -lSDL2


all: $(OBJS)
	$(CC) $(OBJS) $(COMPILE_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)