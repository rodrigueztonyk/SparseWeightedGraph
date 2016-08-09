// SparseWeightedGraph.cpp implements C++ interface for manipulating graphs stored in a sparse format
// Requires nauty to be installed

#include "SparseWeightedGraph.hpp"

// Constructor. Requires number of vertices (_nv) and whether graph is directed (_directed)
SparseWeightedGraph::SparseWeightedGraph(int _nv, bool _directed)
{
	nv = _nv;	// set number of vertices
	nde = 0;	// initialize to 0
	directed = _directed;
	d.resize(nv, 0);	// make d large enough to hold all vertices, initialize to 0
	v.resize(nv, -1);	// make v large enough to hold all vertices, initialize to -1

	if (directed) {
		e.reserve(nv);	// set aside some memory for e
		w.reserve(nv);	// set aside some memory for w
	} else {
		e.reserve(2 * nv);	// set aside some memory for e
		w.reserve(2 * nv);	// set aside some memory for w
	}
}

// Copy constructor
SparseWeightedGraph::SparseWeightedGraph(const SparseWeightedGraph & swg)
{
	nv = swg.nv;
	nde = swg.nde;
	directed = swg.directed;
	nauty_nv = swg.nauty_nv;
	nauty_nde = swg.nauty_nde;
	
	// use std::vector::operator=
	d = swg.d;
	e = swg.e;
	v = swg.v;
	w = swg.w;
	lab = swg.lab;
	ptn = swg.ptn;
	nauty_d = swg.nauty_d;
	nauty_v = swg.nauty_v;
	nauty_e = swg.nauty_e;
	orbits = swg.orbits;
}

// constructor from file. Assumes modified DIMACS format, vertices labeled 0,1,...,(n-1)
SparseWeightedGraph::SparseWeightedGraph(FILE * f)
{
	int c, ch, i, j, m, fileEdges,prev_i=-1,line=-1;
	bool initialized = false;
	char name[40];

	// read in file
	while ((ch = getc(f)) != EOF) {	// while not at the end of the file
		line++;
		switch(ch) {
			case 'p':
				if (!initialized) { // if not initilized yet
					if (fscanf(f, "%s %s %d %d %d", name, name, &nv, &fileEdges, &directed) != 5) { // something ain't right with this file...
						fprintf(stderr, "Something went wrong reading in the graph size. Make sure your first line is right...\n");
						exit(EXIT_FAILURE);
					}
					nde = 0;
					d.resize(nv, 0);	// initialize degrees to 0
					v.resize(nv, -1);	// initialize vs to -1
					if (directed) {
						e.resize(fileEdges,-1);	// make large enough to hold all edges, initialize to -1 to indicate not used yet
						w.resize(fileEdges, 0);	// make large enough to hold all edges, initialize to 0 to indicate not used yet
					} else {
						e.resize(2 * fileEdges, -1);
						w.resize(2 * fileEdges, 0);
					}
					initialized = true;
				} else {
					fprintf(stderr, "We found more than one line beginning with 'p'!\n(This is bad)\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 'c':	// comment line, ignore it
				while ((c = getc(f)) != EOF && c != '\n') continue;
				if (c == EOF) ungetc(c,f);
				break;
			case 'e':	// edge
				if (!initialized) {	// double plus ungood
					fprintf(stderr, "An edge was found before memory was allocated!\nCheck format of input file and ensure it is DIMACS\n");
					exit(EXIT_FAILURE);
				}
				if (fscanf(f, " %d %d %d", &i, &j, &m) != 3) {	// should be 3 items
					fprintf(stderr, "Trouble reading edge on line %d.\nEnsure input file is in DIMACS.\n",line);
					exit(EXIT_FAILURE);
				}
				addEdge(i,j,m);
				if (!directed) addEdge(j,i,m);
				break;
			case '\n':
				break;
			case EOF:
				return;
			default :
				fprintf(stderr, "We found something odd in the file on line %d.\nReally, we've got no idea what to do. Exiting...\n",line);
				exit(EXIT_FAILURE);
		}
	}
}

// Destructor. Should work without doing anything as all objects are std::vector objects
SparseWeightedGraph::~SparseWeightedGraph() { }

// Copy assignment operator
SparseWeightedGraph & SparseWeightedGraph::operator=(const SparseWeightedGraph & swg)
{
	if (this == &swg) return *this;	// if object assigned to itself, do nothing

	// do as in copy constructor
	nv = swg.nv;
	nde = swg.nde;
	directed = swg.directed;
	nauty_nv = swg.nauty_nv;
	nauty_nde = swg.nauty_nde;
	d = swg.d;
	e = swg.e;
	v = swg.v;
	w = swg.w;
	lab = swg.lab;
	ptn = swg.ptn;
	nauty_d = swg.nauty_d;
	nauty_e = swg.nauty_e;
	nauty_v = swg.nauty_v;
	orbits = swg.orbits;
}

// Builds weighted nauty graph
bool SparseWeightedGraph::updateW()
{
	nauty_nv = nv + nde;	// number of nauty vertices
	nauty_nde = 2 * nde;	// number of nauty edges
	nauty_d.resize(nauty_nv,0);	// initialize to 0, unused
	nauty_e.resize(nauty_nde,0);
	nauty_v.resize(nauty_nv,0);
	lab.resize(nauty_nv,0);
	ptn.resize(nauty_nv,0);
	orbits.resize(nauty_nv,0);
	int k = 0;	// used to iterate over edges in weighted graph
	while (k < nde) {
		nauty_e[k] = nv + k;	// bisect arc (i,j)
		nauty_e[nde + k] = e[k];
		k++;
	}

	for (int i = 0; i < lab.size(); i++) {	// begin initial labelling, before sorted by weights
		lab[i] = i;
		orbits[i] = i;
	}

	for (int i = 0; i < nv; i++) {	// set degrees of vertices and places to find adjacencies
		nauty_d[i] = d[i];
		if (v[i] > -1) nauty_v[i] = v[i];
		if (v[i] == -1) nauty_v[i] = 0;
	}

	for (int i = 0; i < nde; i++) {	// set degrees of edge vertices and places to find adj
		nauty_d[nv + i] = 1;
		nauty_v[nv + i] = nde + i;
	}

	std::vector<int> weights = w;	// weight vector to sort
	int w_tp, l_tp;

	// bubble sort to get initial coloring of vertices
	for (int i = 0; i < nde - 1; i++) {
		for (int j = 0; j < nde - 1; j++) {
			if (weights[j] > weights[j + 1]) {
				w_tp = weights[j];
				l_tp = lab[j + nv];
				weights[j] = weights[j + 1];
				lab[j + nv] = lab[j + 1 + nv];
				weights[j + 1] = w_tp;
				lab[j + 1 + nv] = l_tp;
			}
		}
	}

	// give initial coloring, first nv vertices are original vertices, must be in same partition
	for (int i = 0; i < nv - 1; i++) {
		ptn[i] = 1;
	}
	ptn[nv - 1] = 0;

	for (int i = 0; i < nde - 1; i++) {	// determine where new partitions start (different weights)
		if (weights[i + 1] - weights[i]) {
			ptn[nv + i] = 0;
		} else {
			ptn[nv + i] = 1;
		}
	}

	ptn[nauty_nv - 1] = 0;
	return true;
}

// builds the unweighted nauty graph
bool SparseWeightedGraph::updateUW()
{
	nauty_v.resize(nv, 0);
	nauty_e.resize(nde, 0);
	nauty_d.resize(nv, 0);
	lab.resize(nv, 0);
	ptn.resize(nv, 0);
	orbits.resize(nv,0);
	nauty_nv = nv;
	nauty_nde = nde;

	for (int i = 0; i < nv; i++) {
		if (v[i] > -1) nauty_v[i] = v[i];
		if (v[i] == -1) nauty_v[i] = 0;
		nauty_d[i] = d[i];
	}

	for (int i = 0; i < nde; i++) {
		nauty_e[i] = e[i];
	}

	return true;
}

// calls nauty on the given graph
bool SparseWeightedGraph::callNauty(bool _print, bool _trivial)
{
	sparsegraph sg;
	statsblk stats;
	SG_INIT(sg);
	static DEFAULTOPTIONS_SPARSEDIGRAPH(options);
	int min_w = w[0];
	int max_w = w[0];

	for (int i = 0; i < w.size(); i++) {	// find min and max weights
		if (min_w > w[i]) min_w = w[i];
		if (max_w < w[i]) max_w = w[i];
	}

	if (max_w - min_w) {	// if min and max weights are different, true weighted graph
		updateW();
		options.defaultptn = false;
	} else {	// min and max weights the same, can consider as unweighted graph
		updateUW();
	}

	sg.nv = nauty_nv;
	sg.nde = nauty_nde;
	sg.d = &nauty_d[0];
	sg.v = &nauty_v[0];
	sg.e = &nauty_e[0];

	sparsenauty(&sg,&lab[0],&ptn[0],&orbits[0],&options,&stats,NULL);
	// if (_print) printOrbits(_trivial);	// can uncomment when printOrbits is actually added
	return true;
}

// prints all variables
bool SparseWeightedGraph::print()
{
	std::printf("Printing variables!\n");
	std::printf("nv = %d\n",nv);
	std::printf("nde = %d\n",nde);
	std::printf("d = ");
	for (int i = 0; i < d.size(); i++) {
		std::printf("%d ",d[i]);
	}
	std::printf("\ne = ");
	for (int i = 0; i < e.size(); i++) {
		std::printf("%d ",e[i]);
	}
	std::printf("\nw = ");
	for (int i = 0; i < w.size(); i++) {
		std::printf("%d ",w[i]);
	}
	std::printf("\nv = ");
	for (int i = 0; i < v.size(); i++) {
		std::printf("%d ",v[i]);
	}
	std::printf("\nnauty_nv = %d\n",nauty_nv);
	std::printf("nauty_nde = %d\n",nauty_nde);
	std::printf("nauty_d = ");
	for (int i = 0; i < nauty_d.size(); i++) {
		std::printf("%d ",nauty_d[i]);
	}
	std::printf("\nnauty_e = ");
	for (int i = 0; i < nauty_e.size(); i++) {
		std::printf("%d ",nauty_e[i]);
	}
	std::printf("\nnauty_v = ");
	for (int i = 0; i < nauty_v.size(); i++) {
		std::printf("%d ",nauty_v[i]);
	}
	std::printf("\nlab = ");
	for (int i = 0; i < lab.size(); i++) {
		std::printf("%d ",lab[i]);
	}
	std::printf("\nptn = ");
	for (int i = 0; i < ptn.size(); i++) {
		std::printf("%d ",ptn[i]);
	}
	std::printf("\norbits = ");
	for (int i = 0; i < orbits.size(); i++) {
		std::printf("%d ",orbits[i]);
	}
	std::printf("\n");
}

// determine if edge (i,j) is in the graph
bool SparseWeightedGraph::isEdge(int i, int j)
{
	for (int k = v[i]; k < v[i] + d[i]; k++) {
		if (e[k] == j) return true;	// if j is found in i's neighbors
	}
	return false;	// if loop finishes, j isn't a neighbor of i
}

// gets weight of edge (i,j)
int SparseWeightedGraph::getWeight(int i, int j)
{
	for (int k = v[i]; k < v[i] + d[i]; k++) {
		if (e[k] == j) return w[k];	// edge found, return weight
	}
	return 0;	// edge not found, return 0
}

// adds an arc -- no runtime check to see if edge already exists
// !!! USER SHOULD ALWAYS USE addEdge ROUTINE!!!
bool SparseWeightedGraph::addArc(int i, int j, int m)
{
	if (d[i] == 0) {	// we need to find a place to start lit of neighbors
		int k = i;
		// search back until either we hit the beginning or some other list of neighbors
		while (--k >= 0 && d[k] == 0) continue;
		if (k == -1) {	// we hit the beginning without finding another list of neighbors
			v[i] = 0;
		} else {	// we didn't hit the beginning
			v[i] = v[k] + d[k];	// put i's neighbors right after k's to keep things in order
		}

		e.insert(e.begin() + v[i], j);
		w.insert(w.begin() + v[i], m);

		for (int l = i + 1; l < nv; l++) {	// increment the positions of remaining vertices
			if (v[l] > -1) v[l]++;
		}
	} else {	// i already has a list in e
		// keep i's arcs in increasing order, because it is prettier that way
		int k = 0;
		while (k < d[i] && e[v[i] + k] < j) k++;
		e.insert(e.begin() + v[i] + k, j);
		w.insert(w.begin() + v[i] + k, m);
		for (int l = i + 1; l < nv; l++) {	// increment the positions of remaining vertices
			if (v[l] > -1) v[l]++;
		}
	}

	d[i]++;	// one more edge from i
	nde++;	// one more directed edge in graph

	return true;
}

// adds edge (i,j) of weight m to graph. runtime check to make sure edge does not already exist
bool SparseWeightedGraph::addEdge(int i, int j, int m)
{
	if (isEdge(i,j)) {	// this shouldn't happen
		fprintf(stderr,"Tried to add edge (%d,%d), but it is already an edge!\nNothing done.\n",i,j);
		return false;
	} else {	// we're okay to add the edge
		addArc(i,j,m);
	}

	if (!directed) addArc(j,i,m);	// if it isn't a directed graph, need to add other arc too
	return true;
}

// deletes an edge, if it exists
// !!!USER SHOULD ALWAYS USE delEdge ROUTINE!!!
bool SparseWeightedGraph::delArc(int i, int j)
{
	for (int k = 0; k < d[i]; k++) {	// search for j in i's neighbors
		if (e[v[i] + k] == j) {	// j was found, remove it
			e.erase(e.begin() + v[i] + k);
			w.erase(w.begin() + v[i] + k);
			// update where subsequent neighbors start
			for (int l = i + 1; l < nv; l++) {
				if (v[l] > -1) v[l]--;
			}
			d[i]--;	// i has one less neighbor
			if (d[i] == 0) v[i] = -1;	// if i now has no neighbors, make it so
			nde--;
			return true;
		}
	}
	return false;	// edge wasn't found, so wasn't removed
}

// deletes edge (i,j). runtime check to make sure edge is in graph
bool SparseWeightedGraph::delEdge(int i, int j)
{
	if (isEdge(i,j)) {	// this is what we want
		delArc(i,j);
		if (!directed) delArc(j,i);	// if undirected, remove (j,i)
		return true;
	} else {	// ut-oh
		fprintf(stderr, "Tried to remove edge (%d,%d), but it doesn't exist!\nNothing done.\n",i,j);
		return false;
	}
}

// used to change weight of arc (i,j) to m
// !!!USER SHOULD ALWAYS USE changeWeight ROUTINE!!!
bool SparseWeightedGraph::changeArcWeight(int i, int j, int m)
{
	for (int k = v[i]; k < v[i] + d[i]; k++) {
		if (e[k] == j) {	// edge found, change weight
			w[k] = m;
			return true;
		}
	}

	return false;	// edge not found, no change done
}

// used to change weight of arc (i,j) to m. runtime check performed
bool SparseWeightedGraph::changeWeight(int i, int j, int m)
{
	if (changeArcWeight(i,j,m)) {	// try to change arc weight, if successful
		if (!directed && changeArcWeight(j,i,m)) return true;	// if undirected, change weight of (j,i) to m. if successful, return true
	} else {	// tried to change weight, but it didn't work
		fprintf(stderr, "Tried to change weight of edge (%d,%d) to %d, but edge not found!\nNothing done.\n",i,j,m);
		return false;
	}
}

// accessor for number of edges
int SparseWeightedGraph::numEdges()
{
	if(directed) return nde;
	return nde/2;
}

// prints graph in standard format to os
std::ostream & operator<<(std::ostream & os, const SparseWeightedGraph & swg)
{
	os << "p edge directed " << swg.nv << " " << swg.nde << " " << swg.directed << "\n";
	for (int i = 0; i < swg.nv; i++) {
		for (int k = 0; k < swg.d[i]; k++) {
			if ((!swg.directed && i < swg.e[swg.v[i] + k]) || (swg.directed)) {
				os << "e " << i << " " << swg.e[swg.v[i] + k] << " " << swg.w[swg.v[i] + k] << "\n";
			}
		}
	}
	return os;
}
