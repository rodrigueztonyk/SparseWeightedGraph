#include "SparseWeightedGraph.hpp"

int main(int argc, char* argv[])
{
    int i,u,v,t;

    std::ifstream infile(argv[1]);

    std::string line;

    if(infile.good() == false) {
        std::cout<<"Unable to open input file!\n";
        exit(EXIT_FAILURE);
    }

    getline(infile, line);

    std::istringstream in(line);

    int n,m;

    in >> n >> m; // record order and size

    SparseWeightedGraph g(n,false); // create undirected graph g

    for (i = 0; i < m; i++) { // read in edges
        getline(infile, line);
        std::istringstream in(line);
        in >> u >> v >> t;
        t = t + 2; // add 2 to t so it is 2 or 3

        if(g.isEdge(u,v)) { // if edge current exists, don't add it, change weight
            int current_weight = g.getWeight(u,v);
            g.changeWeight(u,v,current_weight*t);
        } else { // if edge doesn't exist, add it
            g.addEdge(u,v,t);
        }
    } // for

    std::cout<<g.numVertices()<<std::endl;
    std::cout<<g.numEdges()<<std::endl;
    return 0;

}
