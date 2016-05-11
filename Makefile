CC=g++

symplectic2: symplectic2.cpp
	$(CC) -DBOOST_LOG_DYN_LINK -Wall -o bin/symplectic2 symplectic2.cpp -lgsl -lgslcblas -lboost_log -lboost_log_setup -lpthread -lboost_thread