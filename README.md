# Triangulator
Implementation of the Bouchitté-Todinca algorithm for finding optimal graph triangulations.

Triangulator supports finding the optimal treewidth and minimum fill-in of a graph, generalized hypertreewidth and fractional hypertreewidth of a hypergraph and total table size of a Bayesian network.


## Use
Build the code with the makefile provided. This should create executable file main.

Triangulator reads the input graph from the standard input and outputs the optimal fill edges or tree decomposition to the standard output. Additional information, like the optimal graph parameter obtained is outputted to the standard error stream.

A basic example would be `./main -minfill < instances/anna.graph`

### Flags
You should use a flag to select the optimized graph measure. Available options are:
* -treewidth : Treewidth of a graph
* -minfill : Minimum fill-in of a graph
* -hyper : Generalized hypertreewidth of a hypergraph
* -frachyper : Fractional hypertreewidth of a hypergraph
* -bayes : Total table size of a Bayesian network

There are also additional flags:
* -pmcprogress : Print information to stderr during the PMC enumeration phase.
* -k : Set an upper bound for the search. Triangulator will report the optimal solution or that the value of the optimal solution is higher than k. Use e.g. `-k=12` Currently works only with treewidth and minfill.
