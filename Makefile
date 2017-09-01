CC=g++
#OPTS=-ggdb
OPTS=-O1
INC=-I./include

#DEF=-DDD_INLINE -Dx86

default: example

example: src/example.o src/su2.o src/Net.o
	$(CC) $(OPTS) src/example.o src/su2.o src/Net.o -o example

src/example.o: src/example.cpp include/Net.h
	echo $(INC)
	$(CC) $(OPTS) $(INC) -c src/example.cpp -o src/example.o 

src/Net.o: src/Net.cpp 
	$(CC) $(OPTS) $(INC) -c src/Net.cpp  -o src/Net.o
	
src/su2.o: src/su2.cpp 
	$(CC) $(OPTS) $(INC) -c src/su2.cpp -o src/su2.o

clean: 
	rm -f src/*.o example
