CC=g++

symplectic: symplectic.cpp
	$(CC) -DBOOST_LOG_DYN_LINK -Wall -o bin/symplectic symplectic.cpp -lgsl -lgslcblas -lboost_log -lboost_log_setup -lpthread -lboost_thread
