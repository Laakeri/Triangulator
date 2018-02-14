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

### File formats
Triangulator supports multiple graph formats. The graphs in folders `instances`, `instances-hyper` and `instances-bayes` give examples of some supported formats.

For normal graphs the formats used in both of the tracks of the [PACE 2017 challenge](https://pacechallenge.wordpress.com/2016/12/01/announcing-pace-2017/) and the dimacs graph format is supported.

For hypergraphs only the format used in the [CSP hypergraph library](https://www.dbai.tuwien.ac.at/proj/hypertree/downloads.html) is supported.

For Bayesian networks the .bif, .dsc and .net formats are supported.
