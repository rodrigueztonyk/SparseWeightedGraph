// SparseWeightedGraph.hpp C++ interface for manipulating weighted graphs stored in a sparse format
// It is assumed vertices are labeled 0,1,...,(n-1)
// This utilizes nauty, so nauty must be installed

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "nausparse.h"
#include "naututil.h"


class SparseWeightedGraph
{
	private:
		bool directed;  // if graph is directed or undirected
		int nv;		// number of vertices
		int nde;	// number of edges
		std::vector<int> d;	// array of degrees. d[i] stores degree of vertex i
		std::vector<int> e;	// list of directed edges
		std::vector<int> w;	// list of weights of edges
		std::vector<int> v;	//v[i] is where vertex i's list of adjacencies start in e
		std::vector<int> lab;	// nauty lab for initial coloring
		std::vector<int> ptn;	// nauty ptn for initial coloring
		int nauty_nv;	// number of vertices in nauty graph
		int nauty_nde;	// number of edges in nauty graph
		std::vector<int> nauty_d;	// array of degrees in nauty graph
		std::vector<int> nauty_e;	// array of directed edges in nauty graph
		std::vector<size_t> nauty_v;	// ith element stores where vertex i's list of adjacencies start in nauty_e
		bool addArc(int i, int j, int m);	// mutator for adding arc (i,j) of weight m
		bool delArc(int i, int j);	// mutator for deleting arc (i,j)
		bool changeArcWeight(int i, int j, int m);	// mutator for changing weight of arc (i,j) to m
		bool updateW();	// builds weighted nauty graph
		bool updateUW();	//builds unweighted nauty graph
	public:
		// constructor that takes number of vertices _nv and whether or not the graph is directed
		// _directed should be true if directed, false if undirected. default value is true
		SparseWeightedGraph(int _nv, bool _directed=true);
		SparseWeightedGraph(const SparseWeightedGraph &);	// copy constructor
		SparseWeightedGraph(FILE * fp);	// constructs graph from file. Must be in modified DIMACS format
		virtual ~SparseWeightedGraph();	// destructor. virtual in case needed as base class
		SparseWeightedGraph & operator=(const SparseWeightedGraph &);	// copy assignment operator
		bool callNauty(bool _print = false, bool _trivial = false);	// used to call nauty. 
		// callNauty has two inputs, _print and _trivial. _print = true will print nauty output to stdout. _trivial = true will print trivial orbits
		bool callNautyForNumber(FILE * _fp);	// used to call nauty to print nauty output to file _fp
		std::vector<int> orbits;	// nauty orbits

		bool print();	// outputs all variables to stdout
		bool printOrbits(bool _print = false);	// prints orbits from nauty to stdout. if _print is true will print all orbits (including trivial), otherwise just prints nontrivial orbits
		int numVertices() const {return nv;}	// accessor for nv
		int numEdges() const {return nde;}	// accessor for nde
		int deg(int i) const {return d[i];}	// accessor for d[i]
		bool isEdge(int i, int j);	// is (i.j) an edge
		int getWeight(int i, int j);	// gets weight of edge (i,j)
		bool addEdge(int i, int j, int m = 1);	// mutator for adding edge (i,j) of weight m. returns true if edge successfully added
		bool delEdge(int i, int j);	// mutator for deleting edge (i,j). returns true if successful
		bool isDirected() const {return directed;}	// accessor for directed
		bool changeWeight(int i, int j, int m);	// mutator for changing weight of edge (i,j) to m. returns true if successful
		int whereInV(int i) const {return v[i];}	// accessor for v[i]

		// access to the beginning of the list of neighbors of i
		std::vector<int>::const_iterator beginNeighbors(int i) const {return (e.begin() + v[i]);}

		// access to an iterator at the end of the list of neighbors of i
		std::vector<int>::const_iterator endNeighbors(int i) const {return (e.begin() + v[i] + d[i]);}

		// access to an iterator of the list of weights of edges from i
		std::vector<int>::const_iterator beginNeighborWeights(int i) const {return (w.begin() + v[i]);}

		//access to an iterator the end of the list of weights of edges from i
		std::vector<int>::const_iterator endNeighborWeights(int i) const {return (w.begin() + v[i] + d[i]);}

		// write graph to stream
		friend std::ostream & operator<<(std::ostream &, const SparseWeightedGraph &);
};
