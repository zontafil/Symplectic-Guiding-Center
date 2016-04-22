CC=g++

symplectic2: symplectic2.cpp
	$(CC) -o bin/symplectic2 symplectic2.cpp systems/*.cpp utils/*.cpp