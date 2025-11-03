#include "headers/ICMole/edge.h"
#include "headers/ICMole/vertex.h"

using namespace std;
using namespace ICMole;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CONSTRUCTORS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
   *  \brief Standard constructor
   *  Can only be called by Graph::addEdge function
   *  \param V1 : First vertex of the edge.
   *  \param V2 : Second vertex of the edge
   *  \param atm : If the edge represent an atom -> address of the atom
   *  \param num : Id of the edge, as defined by the program. Default : 0 (depreciated 0 value)
   *  \param weight : For weighted graph, weight of the edge. Default : 1
   */
Edge:: Edge( Vertex  &nV1,
             Vertex  &nV2,
            const unsigned int& num,
            Graph& par)
            : vertex1(nV1),vertex2(nV2),graph(par),
              bond((Bond*)NULL),weight(1),label(""),num(num)

{

}


/*!
   *  \brief Destructor
   *
   *  Will delete edge pointer from each vertex.
   */
Edge::~Edge() 
{

}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// MISCELLANEOUS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
   *  \fn Vertex&  getOtherVertex(Vertex const &vertex) const
   *  \brief Given a vertex involved in the edge, return the other vertex involved
   *  \param vertex : Vertex to look at
   *  \return Vertex pointer of the other vertex
   *  \throw MoleExcept 1030101 when the given vertex is not part of this edge
   */
 Vertex&  Edge::getOtherVertex(const Vertex  &vertex) const 
{
         if (&vertex == &vertex1) return vertex2;
    else if (&vertex == &vertex2) return vertex1;
    throw MoleExcept(1030201,"Edge::getOtherVertex","Given vertex is not part of this edge");
}


/*!
   *  \fn Vertex& Edge::getOtherVertex(Vertex *const vertex) const
   *  \brief Given a vertex involved in the edge, return the other vertex involved
   *  \param vertex : Vertex to look at
   *  \return Vertex pointer
   *  \throw MoleExcept 1030102 when the given vertex is not part of this edge
   */
 Vertex&  Edge::getOtherVertex(const Vertex *const vertex) const 
{
         if (vertex == &vertex1) return vertex2;
    else if (vertex == &vertex2) return vertex1;
    throw MoleExcept(1030301,"Edge::getOtherVertex","Given vertex is not part of this edge");
}

/*!
   * \fn Vertex*  Edge::shareVertex(Edge const &ed) const
   * \brief Given another edge, search the common vertex and return it
   * \param ed : Edge to look at
   * \return The vertex existing in the two edges or (Vertex*)NULL
   */
const Vertex&   Edge::shareVertex(Edge const &ed) const 
{
  if (&ed.getVertex1() == &vertex1 || &ed.getVertex2()==&vertex1) return vertex1;
  if (&ed.getVertex1() == &vertex2 || &ed.getVertex2()==&vertex2) return vertex2;
   throw MoleExcept(1030401,"Edge::shareVertex","No shared vertex");
}



/*!
   * \fn std::string Edge::toString() const
   *  \brief Returns a description of the edge
   *
   */
std::string Edge::toString() const
{
    ostringstream oss;
    oss << "Num:";
    oss.width(5);
    oss << num;
    oss <<" Ve1: "  << vertex1.getNum()
        <<"<-> Ve2: "<< vertex2.getNum()
        <<"\tW="     << weight
        <<"\tLabel=" << label
        <<"\tParent:"<< &graph
        ;
    return oss.str();
}

/*!
   * \fn void Edge::clone(Edge& ed)
   *  \brief Clone the edge
   *
   *  Copy data of the given edge into this one.
   * \warning It only copy Id, weight and label. It does not copy Vertex pointer nor parent or bond related information.
   * \param ed : Edge to copy data from
   * \param with_parent : copy Bond links. It does not replace the initial edge bond
   */
void Edge::clone(Edge const & ed)
{
    num=ed.num;
    weight=ed.weight;
    label=ed.label;
}
