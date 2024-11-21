

#include <iostream>
//#define ICHEM_DEBUG
#include "headers/ICMole/graphmatch.h"
#include "headers/ICMole/vertex.h"
#include "headers/ICMole/edge.h"

using namespace ICMole;
using namespace std;



GraphMatch::GraphMatch(CliqueData& cliqueData):
    graph(cliqueData.graph),
    cliqueData(cliqueData),
    //matrixbool(0),
    orderON(0),
    orderNO(0),

    minSizeClique(3),
    keepAllCliques(false)
{
graph.renumVertex();

}




/*! \fn void GraphMatch::createMatrix() throw(MoleExcept)
  * Create the adjacency matrix of the graph and order them by decreasing edge size.
  * \exception MoleExcept
  */

void GraphMatch::createMatrix() throw(MoleExcept)
{
    // STEP1 : Check list
    const size_t OrigNVertexs =  graph.numVertex();
    if (OrigNVertexs < minSizeClique)
        throw MoleExcept(1130101,
                         "GraphMatch::createMatrix",
                         "Not enough pairs compared to minimal clique size");

// STEP 2 : Defines original order :
    // orderON helps to make the link between vertex numbering and
    // our filtered numbering
    orderON.clear();
    orderON.reserve(OrigNVertexs);
    // whereas orderNO makes the link between our filtered numbering
    // and the original vertex numbering :
    orderNO.clear();
    orderNO.reserve(OrigNVertexs);


    // posordering is a map associating a number of edges to a vertex
    // The final goal here is to put at the beginning of the matrix
    // vertex that have a maximal number of edges
    // in order to optimize calculation times
    // Beware, its a multimap and not a map, meaning that a same value of
    // number of edges can be associated to one or more vertexs
    std::multimap<size_t,size_t>   posordering;     // [Number of links][pair]
    posordering.clear();

    size_t pos=0;

    for (size_t i=0; i < OrigNVertexs; i++)
    {
        // Given pair has a number of edges lower than the minimal
        // clique size - 1 (-1 is for the given pair)
        // So it cannot be part of the final clique
        const Vertex& vertex=graph.getVertex(i);
#ifdef ICHEM_DEBUG
        cout << vertex.getNum()<< " "<< vertex.numEdges();
#endif
        if (vertex.numEdges() < minSizeClique-1)
        {
#ifdef ICHEM_DEBUG
            cout << "Not enough edges compared to minimal size clique "<<endl;
#endif
            orderON.push_back(-1);
        }
        else
        {
            #ifdef ICHEM_DEBUG
            cout << "pos:"<<pos<<endl;
#endif
            orderON.push_back(pos);
            ++pos;
            // Making link between number of edges and the vertex (original number)
            posordering.insert((pair<size_t,size_t>(vertex.numEdges(),i)));
        }

    }

// STEP 3 - Updating orderNO and orderON by decreasing edge numbers :
    int posNO=0;
    for(multimap<size_t,size_t>::reverse_iterator
        it  = posordering.rbegin();
        it  != posordering.rend();
        it++ ){
        // (*it).second correspond to the original vertex number
        orderNO.push_back((*it).second);
        orderON.at((*it).second)=posNO;
        #ifdef ICHEM_DEBUG
        cout << (*it).first<< " edges for "<< (*it).second<<" => Number "<< posNO<<endl;
#endif
        posNO++;
    }

// STEP 4 : Creating the boolean matrix :
    // Now that we have filtered unwanted vertex and have our final numbering
    // We have the number of vertex that we consider :
    const size_t ReduceNVertexs = orderNO.size();
    // Removing all potential data :
    matrixbool.clear();
    // Reserving memory space for N*N
    matrixbool.reserve(ReduceNVertexs*ReduceNVertexs);
    // And filling them with false values
    matrixbool.assign(ReduceNVertexs*ReduceNVertexs,false);

// STEP 5 : filling in the matrix :
    // Diagonal values (pairs against pairs) are always true :
    for (size_t i=0; i< ReduceNVertexs;++i) matrixbool.at(i*ReduceNVertexs+i)=true;

    // Scanning all edges of the product graph
    for (size_t i=0; i< graph.numEdges();++i)
    {
        const Edge &ed=graph.getEdge(i);

        // Checking the both vertexs are considered here :
        if (orderON.at(ed.getVertex1().getNum())==-1
          ||orderON.at(ed.getVertex2().getNum())==-1)continue;

        // Setting boolean matrix to true :
        matrixbool[orderON.at(ed.getVertex1().getNum())*ReduceNVertexs+orderON.at(ed.getVertex2().getNum())]=true;
        matrixbool[orderON.at(ed.getVertex2().getNum())*ReduceNVertexs+orderON.at(ed.getVertex1().getNum())]=true;
    }


#ifdef ICHEM_DEBUG
    for (size_t i=0; i< ReduceNVertexs;++i)
    {
        cout << i <<"\t "<< orderNO.at(i)<<"\t[";
        for (size_t j=0; j< ReduceNVertexs;++j) cout<<" " <<matrixbool[i*ReduceNVertexs+j];
        cout << " ]"<<endl;
    }
#endif

}



void GraphMatch::runSearchAllMaxCliques(const bool verbose) throw(MoleExcept)
{
    // Number of vertex to consider in the product graph:
    const size_t ReduceNVertexs= orderNO.size();
    if (ReduceNVertexs == 0)
        throw MoleExcept(10000,
                         "GraphMatch::runSearchAllMaxCliques",
                         "Not enough vertex in product graph");
    // List of already considered vertex :
    vector<size_t> alreadyDone; alreadyDone.reserve(ReduceNVertexs);

    vector<bool> Candidates;

    std::vector<int> clique;

    bool allowed=false;

    // minus 1 since there is no need to consider the last vertex
    for (size_t i=0; i<ReduceNVertexs-1;++i)
    {
        // Checking is the given vertex is not part of a former detected clique
       if (std::find(alreadyDone.begin(),
                     alreadyDone.end(),
                     i )
                  != alreadyDone.end()) continue;
#ifdef ICHEM_DEBUG

       cout << "#########################"<<endl
          << "#########################"<<endl
          << "NEW RUN : " << i <<endl;
#endif
       const unsigned int &numRef= orderNO.at(i);
       Candidates.clear();
       allowed=false;

       for (size_t j=0; j < ReduceNVertexs;++j)
       {

          if (j >i && matrixbool[i*ReduceNVertexs+j])
           {
               Candidates.push_back(1);
               alreadyDone.push_back(j);
               allowed=true;
           }
           else Candidates.push_back(0);

       }
        if (!allowed)continue;
       clique.clear();clique.push_back(numRef);
       cliqueSearch(Candidates,clique, 1);
    }

    const vector<vector<int> > &listClique=cliqueData.cliques;
    if (listClique.size()==0 || !verbose)return;

    map<unsigned short,unsigned short> sizes;
    for (vector<vector<int> >::const_iterator it = listClique.begin(); it !=listClique.end();++it)
    {
        const vector<int>&clique_out = *it;
        cout << "CLIQUES : ";
        for (vector<int>::const_iterator itC= clique_out.begin(); itC != clique_out.end() ;++itC)
        {
            cout << (*itC)<<" ";
        }

        cout <<endl;
        if (sizes.find(clique_out.size()) == sizes.end())
        {
            sizes.insert(pair<unsigned short,unsigned short>(clique_out.size(),1));
        }
        else sizes.at(clique_out.size())++;
    }

    ostringstream
            oss1,oss2;
    oss1<<"  ->      CLIQUE SIZE : ";
    oss2<<"  ->         N CLIQUE : ";
    for (map<unsigned short,unsigned short>::iterator it = sizes.begin(); it != sizes.end(); ++it)
    {
        oss1<<(*it).first<<"\t";
        oss2<<(*it).second<<"\t";
    }
    cout << oss1.str()<<endl;
    cout << oss2.str()<<endl;



}



void GraphMatch::cliqueSearch(vector<bool>& Candidates,
                              vector<int> &Clique,
                              const int& level)
{
    // Number of vertex to consider in the product graph:
    const size_t ReduceNVertexs= orderNO.size();

    vector<int> new_clique;
    vector<bool> newCandidates;
    vector<unsigned int> tried;
    bool possible;
unsigned short n_cand=0;

#ifdef ICHEM_DEBUG
    ostringstream oss;
    for ( int i=0; i<level; i++) oss <<"|  "; oss <<"-> ";
    cout << oss.str()<<"   – CANDIDATES : ";
    for (size_t i=0; i< Candidates.size();++i)  cout << Candidates.at(i);
    cout <<endl;

#endif


    for (size_t pos =0; pos < ReduceNVertexs;++pos)
    {
        // If the given vertex is not a candidate => continue
        if (!Candidates.at(pos))continue;


        // If a vertex A is linked to an already checked vertex B, it means
        // that vertex A is already in the clique(s) detected from vertex B
        // So we don't need to try them
        if (!tried.empty())
        {

         possible=true;
         for (vector<unsigned int>::iterator j = tried.begin(); j != tried.end(); j++)
         {

             if (matrixbool[pos*ReduceNVertexs+*j] == true) {possible=false;break;}
         }
         if (!possible)continue;
         for (vector< int>::iterator j = Clique.begin(); j != Clique.end(); j++)
         {

             if (matrixbool[pos*ReduceNVertexs+orderON.at(*j)] == true) {possible=false;break;}
         }
         if (!possible)continue;
        }
        tried.push_back(pos);
        new_clique.clear();
        new_clique=Clique;
        new_clique.push_back(orderNO.at(pos));
#ifdef ICHEM_DEBUG
        cout << oss.str()<<"   – LOOKING AT " << pos <<endl;
#endif
        newCandidates.clear();
        n_cand=0;
        for (unsigned int pos2=0; pos2<=pos;pos2++) newCandidates.push_back(false);
        for (unsigned int pos2=pos+1; pos2<ReduceNVertexs;pos2++)
        {
            if (Candidates.at(pos2) &&
                     matrixbool[pos*ReduceNVertexs+pos2])
            {
                n_cand++;
                newCandidates.push_back(true);
            }else newCandidates.push_back(false);
        }
#ifdef ICHEM_DEBUG
        cout << oss.str()<<"   – Number of candidates :" << n_cand <<endl;
#endif

        if (n_cand==0)
        {
#ifdef ICHEM_DEBUG
 cout <<oss.str()<< "   - END CLIQUE"<<endl;
#endif
            if (new_clique.size() < minSizeClique) continue;
            sort(new_clique.begin(), new_clique.end());
            cliqueData.cliques.push_back(new_clique);

            return;
        }
        else if (keepAllCliques)
        {
            sort(new_clique.begin(), new_clique.end());
            cliqueData.cliques.push_back(new_clique);
            cliqueSearch(newCandidates,new_clique,level+1);
        }
        else
        {
            cliqueSearch(newCandidates,new_clique,level+1);
        }

    }



}


