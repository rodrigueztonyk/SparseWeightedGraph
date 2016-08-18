#include "SparseWeightedGraph.hpp"

void readGraph(char* graph_file, SparseWeightedGraph & swg)
{
    int i,u,v,t,n,m;

    std::ifstream infile(graph_file);

    std::string line;

    if(infile.good() == false) {
        std::cout<<"Unable to open input file\n";
        exit(EXIT_FAILURE);
    }

    getline(infile,line);
    std::istringstream in(line);

    in >> n >> m;

    SparseWeightedGraph g(n,false);

    for (i = 0; i < m; i++) {
        getline(infile,line);
        std::istringstream in(line);
        in >> u >> v >> t;
        t += 2;

        if (g.isEdge(u,v)) { // don't add edge, update weight
            g.changeWeight(u,v,g.getWeight(u,v) * t);
        } else { // add edge
            g.addEdge(u,v,t);
        }
    }
    swg = g;
}

int main(int argc, char* argv[])
{

    SparseWeightedGraph g(1);
    readGraph(argv[1],g);

    g.callNauty(false,false);

    std::cout<<"Group size : "<<g.groupSize()<<std::endl;

    return 0;

}
