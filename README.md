# Triangulator
Implementation of the Bouchitté-Todinca algorithm for finding optimal graph triangulations.

Triangulator supports finding the optimal treewidth and minimum fill-in of a graph, generalized hypertreewidth and fractional hypertreewidth of a hypergraph and total table size of a Bayesian network.


## Use
Build the code with the makefile provided. This should create executable file main.

Triangulator reads the input graph from the standard input and outputs the optimal fill edges or tree decomposition to the standard output. Additional information, like the optimal graph parameter obtained is outputted to the standard error stream.

A basic example would be `./main -minfill < instances/anna.graph`
