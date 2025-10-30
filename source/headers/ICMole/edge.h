#ifndef EDGE_H
#define EDGE_H

#include "global.h"

namespace ICMole
{


class Edge
{

    friend class Graph;
private:
    Vertex& vertex1;        /*!< First vertex involved in the edge*/
    Vertex& vertex2;        /*!< Second vertex involved in the edge*/
    Graph& graph;         /*!< Parent of the edge*/
     Bond* bond;           /*!< Possibly related bond */

        double  weight;         /*!< Weight of the edge, if necessary. Default 1 */
   std::string  label;          /*!< Label of the edge*/
   unsigned int  num;            /*!< Id of the edge, as defined by the program*/




////// Constructors :
           Edge(Edge const &);
           Edge& operator=(Edge const &);/* Copy constructor private -> Unavailable */
           Edge( Vertex  &V1,   Vertex  &V2,  const unsigned int& num,Graph&  par);

           ////// Destructors :
                       ~Edge() ;
public:



///// Getters:
  unsigned int  getNum()     const {return num;    }   /*!< \brief Return the Id of the edge */
        double  getWeight()  const {return weight; }   /*!< \brief Return the weight of the edge */
 Vertex& getVertex1() const {return vertex1;}   /*!< \brief Return the first vertex involved in the bond */
  Vertex& getVertex2() const {return vertex2;}   /*!< \brief Return the second vertex involved in the bond */
         Graph& getGraph()   const {return graph; }   /*!< \brief Return the parent graph */
   std::string  getLabel()   const {return label;  }   /*!< \brief Return the label of the edge */


///// Setters
          void  setNum   (const unsigned int& n)  {   num = n; }        /*!< \brief Set the Id of the edge (as given by the Graph) */
          void  setWeight(const       double& n)  {weight = n; }        /*!< \brief Set the weight of the edge */
          void  setLabel (const std::string &lab) { label = lab; }        /*!< \brief Set the label of the graph */

///// Miscellaneous
         Vertex &getOtherVertex(Vertex  const &vertex) const ;
         Vertex &getOtherVertex(const Vertex * const vertex) const ;
        const Vertex &shareVertex(Edge    const &edge  ) const  ;
   std::string  toString() const;
          void  clone(Edge const & ed);
          void setBond( Bond* newbond){bond=newbond;}
          Bond* getBond()  {return bond;}
};


}
#endif // EDGE_H
