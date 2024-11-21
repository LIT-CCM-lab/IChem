#ifndef GRAPH_H
#define GRAPH_H

#include "global.h"

namespace ICMole
{


class Graph
{
private:
    Graph(Graph const &);           // No definition
    Graph& operator=(Graph const &);// No definition



protected:
    VertexList  vertexs;          /*!< \brief List of vertexs - Memory automatically handled */
      EdgeList  edges;            /*!< \brief List of edges - Memory automatically handled */
  unsigned int  maxNumVe;         /*!< \brief Maximum Id for vertex - Automatically updated when adding or deleting a vertex */
  unsigned int  maxNumEd;         /*!< \brief Maximum Id for edge - Automatically updated when adding or deleting an edge */

public:

  ////// Constructors :
      Graph(const unsigned int &NVertex=50, const unsigned int &NEdges=100);

  ////// Destructors :
      ~Graph();
  ////// Vertex :
      Vertex& addVertex(const double &weight=1, const std::string &label="") throw(MoleExcept);
        void  addVertexs(const unsigned int &NtoAdd)                         throw(MoleExcept);
        void  delVertex( Vertex *const vertex, const bool& with_maxnum=true)  throw(MoleExcept);
        void  delVertexs( VertexList& vertexs, const bool& renum=true)       throw(MoleExcept);
        void  renumVertex();
        void  reserveVertex(const size_t& size){vertexs.reserve(size);}
        void  reserveEdges(const size_t& size){edges.reserve(size);}
const Vertex& getVertex(const size_t &n,const bool &pos=true) const    throw(MoleExcept);
unsigned int  getMaxNumVert()   const  {return maxNumVe;}                                       /*!< \brief Return the maximal Num for Vertexs */
size_t  numVertex()       const  {return vertexs.size();}                   /*!< \brief Return the number of vertex in this graph */
        ItVert firstVertex() {return vertexs.begin();}
        ItVert lastVertex() {return vertexs.end();}
        void clear();

  ////// Edge :
        Edge& addEdge( Vertex &vertex1, Vertex &vertex2) throw(MoleExcept);
        void  delEdge(const Edge* const   edge)          throw(MoleExcept);
        void  delEdge(const   EdgeList& EdList)          throw(MoleExcept);
        void  delEdge(const       Edge &  edge)          throw(MoleExcept);
        void  renumEdge();
unsigned int  getMaxNumEdge()   const  {return maxNumEd;}                       /*!< \brief Return the maximal Num for edges */
const Edge& getEdge(const size_t &pos) const throw(MoleExcept);
size_t numEdges()        const  {return edges.size();}     /*!< \brief Return the number of edges in this graph */
ItEdge firstEdge() {return edges.begin();}
ItEdge lastEdge() {return edges.end();}

        void  createCopyFrom(const Graph &Gr);

  ///// Miscellaneous:
 std::string  toString() const;
      double  getDensity() const ;
      void getDistMatrix(std::vector<unsigned short> &matrix);

};


}
#endif // GRAPH_H
