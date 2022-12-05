CC=g++ -march=native -O3
CFLAGS=-c -I. -std=c++0x -Wfatal-errors

all: TC_index TC_query

TC_index: TreeCount.o Graph.o TreeDecomp.o
	$(CC) TreeCount.o Graph.o TreeDecomp.o -o TC_index
	
TC_query: TreeCountQuery.o Graph.o TreeDecomp.o
	$(CC) TreeCountQuery.o Graph.o TreeDecomp.o -o TC_query

TreeCountQuery.o: TreeCountQuery.cpp
	$(CC) $(CFLAGS) TreeCountQuery.cpp -o TreeCountQuery.o

TreeCount.o: TreeCount.cpp
	$(CC) $(CFLAGS) TreeCount.cpp -o TreeCount.o

Graph.o: Graph.cpp
	$(CC) $(CFLAGS) Graph.cpp -o Graph.o

TreeDecomp.o: TreeDecomp.cpp
	$(CC) $(CFLAGS) TreeDecomp.cpp -o TreeDecomp.o

.PHONY: clean
clean:
	-rm *.o
