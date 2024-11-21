#ifndef VERTEX_H
#define VERTEX_H

#include "global.h"


namespace ICMole
{

class Vertex
{
    friend class Edge;
    friend class Graph;

    unsigned int  num;       /*!< \brief  Number as assigned by the program */
          double  weight;    /*!< \brief  Weight for the vertex*/
     std::string  label;     /*!< \brief  Label of the vertex */
           Graph* graph;    /*!< \brief  Parent graph*/
            Atom* atom;
        EdgeList  links;     /*!< \brief  List of edges */
        Vertex(Vertex const &);// No definition
        Vertex& operator=(Vertex const &);// No definition

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ///////////////// CONSTRUCTORS /////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////

     Vertex(Graph * const parent,
            const unsigned int &num,
            const double &weight=1,
            const std::string &label ="");

     ~Vertex();

     ////////////////////////////////////////////////
     ////////////////////////////////////////////////
     ///////////////// PRIVATE EDGE /////////////////
     ////////////////////////////////////////////////
     ////////////////////////////////////////////////

             void  addEdge(Edge *const ed);
             void  delEdge(const Edge *const ed) throw(MoleExcept);



public :

///// Edges :
void reserve(const size_t& size){links.reserve(size);}
             bool  hasEdgeWith(Vertex const &ve) const;//
       const Edge* getEdgeWith(Vertex const &ve) const throw(MoleExcept);//
             void  cleanEdge()                         throw(MoleExcept);
        size_t numEdges() const {return links.size();}
       const Edge* getEdge(const size_t& i) const throw(MoleExcept);
         Vertex& getVertex(const size_t& i) const throw(MoleExcept);

///// Setters:
             void  setParent(      Graph  *const  gr  ){  graph=gr;  }       /*!< \brief Set the parent graph */
             void  setNum   (const unsigned int &nnum=0){    num=nnum;}       /*!< \brief Set the Id of the vertex */
             void  setWeight(const double       &wei =0){ weight=wei; }       /*!< \brief Set the Weight of the vertex*/
             void  setLabel (const std::string  &str   ){  label=str; }       /*!< \brief Set the label of the vertex*/

///// Getters:
const unsigned int&  getNum()    const { return num;   }/*!< \brief Return the id of the vertex, as given by the program*/
           double  getWeight() const { return weight;}/*!< \brief Return the weight of the vertex*/
      std::string  getLabel()  const { return label; }/*!< \brief Return the label of the vertex*/
            Graph& getParent() const { return *graph;}

//// Miscellaneous :
             void  clone(const Vertex& vertex);
      std::string  toString(const bool &no_edge=true) const;
             void  setAtom( Atom* const atm) {atom=atm;}
             Atom* getAtom()const {return atom;}
};

}


#endif // VERTEX_H
