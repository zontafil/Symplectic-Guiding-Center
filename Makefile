CC=g++

symplectic2: symplectic2.cpp
	$(CC) -Wall -o bin/symplectic2 symplectic2.cpp -lgsl -lgslcblas