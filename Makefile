CC=g++

symplectic: symplectic.cpp
	$(CC) -g -DBOOST_LOG_DYN_LINK -Wall -o bin/symplectic symplectic.cpp emfields/eqdskReader/eqdskPsiInterp.cpp emfields/eqdskReader/eqdskReader.cpp emfields/ascot5-spline/interp2Dexpl.c emfields/ascot5-spline/interp1Dcomp.c emfields/ascot5-spline/spline1Dcomp.c emfields/ascot5-spline/splineexpl.c -lgsl -lgslcblas -lboost_log -lboost_log_setup -lpthread -lboost_thread
