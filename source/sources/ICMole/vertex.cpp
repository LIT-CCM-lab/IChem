#include <iterator>
#include <algorithm>
#include "headers/ICMole/vertex.h"
#include "headers/ICMole/edge.h"
#include "headers/ICMole/graph.h"
#include "headers/ICMole/atom.h"

/*!
\todo Add (Atom*)NULL as initializer of Vertex
*/

using namespace ICMole;
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// CONSTRUCTORS /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*! \fn Vertex::Vertex(  Graph * const parent, const unsigned int num,const double weight, const std::string label)
  * \brief Standard constructor
  * \param parent : Parent graph related to this vertex
  * \param num    : Id of the vertex (as given by the program)
  * \param weight : Weight of the vertex: Default: 0
  * \param label  : Label of the vertex: Default: empty
  * \warning : Private function, can only be called by either Graph, Atom or molecule object
  *
  *
  */
Vertex::Vertex(   Graph      * const parent,
                 const unsigned int &num,
                 const       double &weight,
                 const  std::string &label):
    num(num),weight(weight),label(label),graph(parent),atom((Atom*)NULL)
{
links.reserve(4);
}


Vertex::~Vertex()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////// EDGES /////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*! \fn void  Vertex::addEdge     (const Edge  * const ed)
  * \brief Adding an edge to the vertex
  * \param ed : Edge to add
  */
void  Vertex::addEdge    ( Edge  *const  ed)
{

      links.push_back(ed);
}


/*! \fn void  Vertex::delEdge    (const  Edge  *const ed )
  * \brief Delete an edge to the vertex edge list
  * \param ed : Edge to remove
  * Check whether the edge is in the edge list of the vertex and delete it if found<br/>
  * \warning It doesn't delete the edge, just the link in the vertex.
  * \throw MoleExcept 1020101 when the given edge is not part of this vertex
  */
void  Vertex::delEdge    (const  Edge  *const ed ) throw(MoleExcept)
{
    const ItEdge ited = find(links.begin(),links.end(),ed);
    if (ited == links.end())
        throw MoleExcept(1020101,
                         "Vertex::delEdge",
                         "Given edge is not part of this vertex : \n"+toString()+"\n"+ed->toString());
    links.erase(ited);
}



/*! \fn bool  Vertex::hasEdgeWith(Vertex const &ve)
  * \brief Check if the given vertex shares an edge with this vertex
  * \param ve : Vertex to check
  * \return true when an edge exists. False otherwise
  */
bool  Vertex::hasEdgeWith(Vertex const &ve) const
{
    const Edge * ed=(Edge*)NULL;
    for (ItCEdge it = links.begin(); it != links.end(); it++)
      {
        ed = *it;
        if (&ed->getVertex1() == &ve || &ed->getVertex2()==&ve) return true;
      }
    return false;
}




/*! \fn const Edge* Vertex::getEdgeWith(Vertex const &ve)
  * \brief Check if the given vertex shares an edge with this vertex, and returns the corresponding edge
  * \param ve : Vertex to check
  * \return Pointer of the edge betwwen ve and this vertex, or (Edge*)NULL when no edge exists.
  * \throw 1030301 - getOtherVertex
  */
const Edge* Vertex::getEdgeWith(Vertex const &ve) const throw(MoleExcept)
    {
    try
    {
  const Edge *ed;
  for (ItCEdge it = links.begin(); it != links.end(); it++)
    {
      ed = *it;
      if (&(ed->getOtherVertex(this)) == &ve) return ed;
    }
        return (Edge*)NULL;
    }catch (MoleExcept &e)
    {
        //1030301 - getOtherVertex
        e.addTrace("Vertex::getEdgeWith");
        throw;
    }
    }


/*! \fn void Vertex::cleanEdge() throw(MoleExcept)
  * \brief Delete all edges of this vertex
  * Calls delEdge() function of the parent graph to make a clean deletion of all edges made with this vertex
  * \throw See Graph::delEdge(). If 1040201: Fatal error
  */
void Vertex::cleanEdge() throw(MoleExcept)
{

    try
    {
        graph->delEdge(links);
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Vertex::cleanEdge()\nFrom : "+toString());
    }

      links.clear();

}


const Edge* Vertex::getEdge(const size_t& i) const throw(MoleExcept)
{
    if (i >= links.size()) throw MoleExcept(1020201,
                                                   "Vertex::getEdge",
                                                   "Value above the number of edges");
    return links.at(i);
}

 Vertex& Vertex::getVertex(const size_t& i) const throw(MoleExcept)
{
    if (i >= links.size()) throw MoleExcept(1020201,
                                                   "Vertex::getEdge",
                                                   "Value above the number of edges");
    return links.at(i)->getOtherVertex(this);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////// MISCELLANEOUS /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \fn std::string Vertex::toString(const bool no_edge) const
  * \brief Output verbose information about this vertex
  * \param no_edge : Default TRUE. Also output linked edge information
  */
std::string Vertex::toString(const bool &no_edge) const
{
    ostringstream oss;


    oss << "VERTEX - ";
    oss.width(5);oss<< num <<". W="<< weight;
    if (label.length())         oss <<". Label:"<< label;
    if (atom != (Atom*)NULL) oss << ". Atom:"<<atom->getIdentifier();

    if (links.size() && !no_edge){
            oss <<endl;
        try
        {
        for (ItCEdge it = links.begin(); it != links.end(); it++)
          {
            const Vertex& vert2=(*it)->getOtherVertex(*this);
            oss << "   |->Edge with ("<<vert2.getNum()<<". W="<<vert2.getWeight()<<") W="<<(*it)->getWeight()<<endl;
          }
        }
        catch (MoleExcept &e)
        {
            oss << "Issue while searching edge information"<<endl<< e.getCode()<<endl<<e.getData()<<endl<<e.getSource()<<endl<<e.getTrace()<<endl;
        }
       }

    return oss.str();
}



/*! \fn void Vertex::Clone(Vertex& vertex)
  * \brief Clone the data of the given vertex to this vertex
  * \param vertex : Vertex to copy data from
  * \param with_linker : Default false. If true, also copy adress of related atom;
  * Clone Ids, weight and labels of the given vertex to the current one. However, it does not copy  links.
  * Therefore no edge will be created.
  */
void Vertex::clone(const Vertex& vertex)
    {
           num = vertex.num;
        weight = vertex.weight;
         label = vertex.label;
    }

