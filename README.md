# Triangulator
Implementation of the Bouchitté-Todinca algorithm for finding optimal graph triangulations.

Triangulator can be used for finding an optimal solution for the following problems: The treewidth and minimum fill-in of a graph, the generalized hypertreewidth and fractional hypertreewidth of a hypergraph, the total table size of a Bayesian network.


## Use
Build the code with the makefile provided. This should create executable file main.

Triangulator reads the input graph from the standard input and outputs the optimal fill edges or tree decomposition to the standard output. Additional information, like the optimal graph parameter obtained is outputted to the standard error stream.

A basic example would be `./main -minfill < instances/anna.graph`

### Flags
You should use a flag to select the graph parameter that is optimized. Available options are:
* -treewidth : Treewidth of a graph
* -minfill : Minimum fill-in of a graph
* -hyper : Generalized hypertreewidth of a hypergraph
* -frachyper : Fractional hypertreewidth of a hypergraph
* -bayes : Total table size of a Bayesian network

There are also additional flags:
* -pmcprogress : Print information to stderr during the PMC enumeration phase.
* -k : Set an upper bound for the search. Triangulator will report the optimal solution or that the value of the optimal solution is higher than k. Use e.g. `-k=12` Currently works only with treewidth and minimum fill-in.

### File formats
Triangulator supports multiple graph formats. The graphs in directories `instances`, `instances-hyper` and `instances-bayes` give examples of some supported formats.

For normal graphs the formats used in both of the tracks of the [PACE 2017 challenge](https://pacechallenge.wordpress.com/2016/12/01/announcing-pace-2017/) and the [DIMACS](http://prolland.free.fr/works/research/dsat/dimacs.html) graph format are supported.

For hypergraphs the format used in the [CSP hypergraph library](https://www.dbai.tuwien.ac.at/proj/hypertree/downloads.html) is supported.

For Bayesian networks the .bif, .dsc and .net formats are supported.

### Additional notes
To use GLPK as the LP solver for Fractional hypertreewidth you need to modify `src/Makefile`.

The output formats for hypergraphs and Bayesian networks are likely not very useful. Feel free to open an issue if there is an output format that would be useful for your project.
