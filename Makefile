CC = g++
CFLAGS = -std=c++11 -O2 -Wall -Wextra -Wshadow -g -Wfatal-errors

solver:
	$(MAKE) -C src all
	cp src/solver main

clean: 
	$(MAKE) -C src clean
	rm -f main
