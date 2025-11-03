#ifndef GRAPHMATCH_H
#define GRAPHMATCH_H

#include "headers/ICMole/graph.h"

namespace ICMole
{

struct Pair
{
    const  Atom*  ref;
    const  Atom*  comp;
    Vertex*   vertex;
    double weight;
    unsigned short PType;
    Pair(const Atom* ref,
         const Atom* comp,
         Vertex* vertex,
         const double& weight=-1,
         const unsigned short& PType=0):
        ref(ref),
        comp(comp),
        vertex(vertex),
        weight(weight),
        PType(PType){}
    Pair(const Pair& p):
        ref(p.ref),comp(p.comp),vertex(p.vertex),weight(p.weight),PType(p.PType)
    {

    }
};

struct CliqueData
{
    std::vector<Pair> listPairs;
    Graph graph;
    std::vector< std::vector<int> > cliques;
};


class GraphMatch
{
private:
    /**
     * @brief graph : Graph product used for clique detection
     */
    Graph &graph;


    /**
     * @brief cliqueData : data structure containing the graph,
     * the list of pairs and the clique results
     */
    CliqueData & cliqueData;

    /**
     * @brief matrixbool : boolean matrix for edge linking
     */
    std::vector<bool> matrixbool;

    /**
     * @brief orderON : from original graph to newer graph, makes
     * the mapping of vertex ids
     */
    std::vector<int> orderON;

    /**
     * @brief orderNO : from newer graph to older graph, makes
     * the mapping of vertex ids
     */
    std::vector<int> orderNO;


    /**
     * @brief minSizeClique : minimal size of the clique
     */
    size_t minSizeClique;

    /**
     * @brief keepAllCliques
     */
    bool keepAllCliques;

public:
    GraphMatch(CliqueData&);
    void createMatrix() ;
    void runSearchAllMaxCliques(const bool verbose) ;
    void cliqueSearch(std::vector<bool>& Candidates,
                      std::vector<int> & Clique,
                      const int& level);
    void setMinSizeClique(const size_t& newSize){ minSizeClique=newSize;}
};

}
#endif // GRAPHMATCH_H
