LPSOLVER = ADHOC

CC = g++
CFLAGS = -std=c++11 -O2 -Wall -Wextra -Wshadow -g -Wfatal-errors
OBJFILES = bayes.o graph.o hyper.o IO.o MCS.o minfill.o PMC.o setcover.o timer.o treewidth.o
LIBS = 

ifeq ($(LPSOLVER), GLPK)
	CFLAGS += -DUSE_GLPK -I../glpk-4.63/incl -L../glpk-4.63/incl
	LIBS += -lglpk
endif

all: solver

solver: solver.o $(OBJFILES)
	$(CC) $(CFLAGS) solver.o $(OBJFILES) -o solver $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -MMD -c $< -o $@ $(LIBS)

-include *.d

clean: 
	rm -f *.o *.d