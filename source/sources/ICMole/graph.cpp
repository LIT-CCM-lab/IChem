#include <iostream>
#include <stdexcept>
#include "headers/ICMole/graph.h"
#include "headers/ICMole/vertex.h"
#include "headers/ICMole/edge.h"
#include "headers/ICMole/atom.h"


using namespace std;
using namespace ICMole;

/*! \fn Graph::Graph(const unsigned int NVertex, const unsigned int NEdges)
  * \brief Constructor
  * \param NVertex : Number of vertex of this graph. Can be 0. It is used to allocate memory and makes the code faster. Default:50
  * \param NEdges : Number of edges of this graph. Can be 0. It is used to allocate memory and makes the code faster. Default:100
  */

Graph::Graph(const unsigned int &NVertex, const unsigned int &NEdges):maxNumVe(0),maxNumEd(0)
{
    vertexs.reserve(NVertex);
    edges.reserve(NEdges);
}


Graph::~Graph()
{
    for(size_t i = 0; i < edges.size(); ++i)   delete edges[i];
    for(size_t i = 0; i < vertexs.size(); ++i) delete vertexs[i];
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// VERTEXS ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
  * \fn Vertex&  Graph::addVertex(const double weight, const std::string label) 
  * \param weight: Weight of the vertex- Default 1
  * \param label : Label of the vertex- Default empty
  * \brief Add the given vertex to the graph
  * \return The newly created vertex
  * \throw 1040101 Bad allocation exception
  */
Vertex&  Graph::addVertex(const double &weight, const string &label) 
{
    Vertex *vertex= (Vertex*)NULL;
    try{
        vertex= new Vertex(this,maxNumVe,weight,label);
        vertexs.push_back(vertex);
        maxNumVe++;
        return *vertex;
    }
    catch (std::bad_alloc& ba)
    {
        if (vertex != (Vertex*)NULL)delete vertex;
        string error("Bad allocation exception - exiting : ");
        error+=ba.what();
        throw MoleExcept(1040101,"Graph::addVertex",error);
    }
}

/*! \fn void Graph::addVertexs(const unsigned int N)
 * \brief Creates N vertexs within this Graph
 * \param N : Number of vertex to be created
 * \throw 1040102 : When bad allocation appends
 */

void Graph::addVertexs(const unsigned int &N) 
{
    if (N==0)return;
    if (N+vertexs.size() > vertexs.capacity()) vertexs.reserve(N+vertexs.size());
    Vertex *vertex= (Vertex*)NULL;
    try
    {
    for (unsigned int i=0; i<N;i++)
    {
        vertex= new Vertex(this,maxNumVe);
        vertexs.push_back(vertex);
        maxNumVe++;
        vertex=(Vertex*)NULL;
    }
    }
    catch (std::bad_alloc& ba)
    {
        if (vertex != (Vertex*)NULL)delete vertex;
        string error("Bad allocation exception - exiting : ");
        error+=ba.what();
        throw MoleExcept(1040102,"Graph::addVertexs",error);
    }
}



/*! \fn void Graph::delVertex(Vertex* const vertex, const bool with_maxnum)
  * \brief Remove the given vertex from the graph.
  * It will first search in the graph if this vertex is associated to this graph, then delete the vertex
  * It also call the cleanEdge() function of the vertex and therefore delete all the associated edges of this graph
  * Also update the maxNumVe value
  * \throw 1040301 - No vertex given
  * \throw 1040302 - If the graph doesn't have the given vertex
  * \throw Exceptions from Graph::delEdge() && Graph::cleanEdge()
  * \param vertex : Vertex to delete
  * \param with_maxnum : Update the maximum number for vertex (strong recommended)
  * \warning Will call cleanEdge() function even if the given vertex is not in the graph
  */
void Graph::delVertex( Vertex *const vertex, const bool& with_maxnum) 
{
    try{
        if (vertex == (Vertex*)NULL)throw MoleExcept(1040301,"Graph::delVertex","No vertex given");

// Searching the position of the vertex within the vertex list :
        const ItVert itVER=find(vertexs.begin(), vertexs.end(),vertex);
        if (itVER == vertexs.end())throw MoleExcept(1040302,"Graph::delVertex","This graph doesn't contains the following vertex : "+vertex->toString());

// Deleting all edges related :

        vertex->cleanEdge();

        const unsigned int numVe= vertex->getNum();

// Deleting the edge from Graph :
        delete *itVER;
        vertexs.erase(itVER);
        if (!with_maxnum)return;

// Updating the maxNumVe - Unecessary if the vertex num is not the maxNumVe
        if (numVe+1 != maxNumVe) return;
        bool zer=false;
        maxNumVe=1;
        for(ItVert it2 = vertexs.begin(); it2 != vertexs.end(); it2++)
        {
            if ((*it2)->getNum() >= maxNumVe){ maxNumVe= (*it2)->getNum();if (maxNumVe==1) zer=true;}
        }
        if (maxNumVe >1|| zer==true) maxNumVe++;
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Graph::delVertex()");
        switch (e.getCode())
        {
        case 1040701:
        case 1040702:
        case 1040202:
        case 1040101:
            switch (Moleaccess)
            {
            case Levels::NOTICE:
            case Levels::WARNING:cerr <<e.getCode()<<"\t"<<e.getData()<<" in "<<e.getSource()<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:cerr <<e.getCode()<<"\t"<<e.getData()<<" in "<<e.getSource()<<endl; e.setCode(999999);throw;
            }
        }

    }

}


void Graph::delVertexs( VertexList& vertexsToDel, const bool& renum) 
{
    try{

        for (size_t i =0; i < vertexsToDel.size(); ++i)
        {
        Vertex* vertex = vertexsToDel.at(i);
        if (vertex == (Vertex*)NULL)throw MoleExcept(1040401,"Graph::delVertexs","No vertex given");

         ItVert itVER=find(vertexs.begin(), vertexs.end(),vertex); // Deletion if it finds in vertexs the corresponding vertex
        if (itVER == vertexs.end())throw MoleExcept(1040402,"Graph::delVertexs","This graph doesn't contains the following vertex : "+vertex->toString());
vertex->cleanEdge();

        // Deleting the edge from Graph :
        delete *itVER;
        vertexs.erase(itVER);
        }
        if (!renum)return;

        // Updating the maxNumVe - Unecessary if the vertex num is not the maxNumVe
        bool zer=false;
        maxNumVe=1;
        for(ItVert it2 = vertexs.begin(); it2 != vertexs.end(); it2++)
        {
            if ((*it2)->getNum() >= maxNumVe){ maxNumVe= (*it2)->getNum();if (maxNumVe==1) zer=true;}
        }
        if (maxNumVe >1|| zer==true) maxNumVe++;
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Graph::delVertex()");
        switch (e.getCode())
        {
        case 1040701:
        case 1040702:
        case 1040202:
        case 1040101:
            switch (Moleaccess)
            {
            case Levels::NOTICE:
            case Levels::WARNING:cerr <<e.getCode()<<"\t"<<e.getData()<<" in "<<e.getSource()<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:cerr <<e.getCode()<<"\t"<<e.getData()<<" in "<<e.getSource()<<endl; e.setCode(999999);throw;
            }
        }

    }

}




/*! \fn void Graph::renumVertex()
  *\brief Renum all vertexs of the graph
  *
  * Scan all vertex of the graph and renumerotate all. Call setNum() function for all vertexes
  */
void Graph::renumVertex()
{
  maxNumVe=0;
  for (ItVert it = vertexs.begin(); it != vertexs.end(); it++)
    {
      ((*it))->setNum(maxNumVe);
      maxNumVe++;
    }
}


/*! \fn Vertex& Graph::getVertex(const unsigned int n,bool pos)
 *  \brief Return the vertex associated with the number n.
 *  \param pos: When pos=true (default) : n corresponds to the position within the \
 *  list of vertexs. Otherwise if the id thought getNum() function that will be checked to n
 *  \param n : Position
 *  \return The corresponding vertex
 *  \throw 1040801- When pos is true and n is above the number of vertex in the graph
 *
 */
const Vertex& Graph::getVertex(const size_t &n, const bool &pos) const 
{
  if (pos)
    {
      if (n >= vertexs.size()) throw MoleExcept(1040501,"Graph::getVertex"," given value above the number of vertexs within the graph");
      return *vertexs.at(n);
    }
  for (size_t i=0; i<vertexs.size();i++) if (vertexs.at(i)->getNum() == n) return *vertexs.at(i);

  throw MoleExcept(1040502,"Graph::getVertex"," No vertex found with the given id");

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// EDGES ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
  * \fn Edge& Graph::addEdge(const Vertex &ve1,const Vertex &ve2) 
  * \param ve1 : First vertex of the edge
  * \param ve2 : Second vertex of the edge
  * \throw MoleExcept code 1040601 when ve1 and ve2 are the same vertex. Returns the last added edge
  * \throw MoleExcept code 1040602 when ve1 are not part of this graph. Returns the last added edge
  * \throw MoleExcept code 1040603 when ve2 are not part of this graph. Returns the last added edge
  * \throw MoleExcept code 1040604 when memory allocation of the edge failed : Fatal error
  * \return Reference of the newly created edge
  * Add an edge to edges list. <br/>
  * Each time a bond is added to a molecule, the edge associated to the bond is added to the graph associated to the molecule
  */
Edge& Graph::addEdge(Vertex &vertex1, Vertex &vertex2) 
{
    Edge *ed=(Edge*)NULL;
    try
    {
        if (&vertex1 == &vertex2)         throw MoleExcept(1040601, "Graph::addEdge", "Vertexs given for the creation of the edge are the same");
        if (&vertex1.getParent() != this) throw MoleExcept(1040602, "Graph::addEdge", "First Vertex is not part of this graph");
        if (&vertex2.getParent() != this) throw MoleExcept(1040603, "Graph::addEdge", "Second Vertex is not part of this graph");
        ed = new Edge(vertex1,vertex2,maxNumEd,*this);
        vertex1.links.push_back(ed);
        vertex2.links.push_back(ed);
        maxNumEd++;
        edges.push_back(ed);
        return *ed;
    }
    catch(MoleExcept &e)
    {
        if (ed!=(Edge*)NULL) delete ed;
        throw;
    }
    catch (std::bad_alloc& ba)
    {
        Moleaccess=Levels::FATAL;
        string error("Bad allocation exception - exiting :"); error+=ba.what();
        throw MoleExcept(1040604,"Graph::addEdge", error);
    }
}






void Graph::delEdge(const Edge* const edge) 
{
    if (edge == (Edge*)NULL) throw MoleExcept(1040701,"Graph::delEdge","No edge given");
    const unsigned int NumEd= edge->getNum();
// Check edge existence
    const ItEdge itPED=find(edges.begin(), edges.end(),edge); // Deletion if it finds in edges the corresponding edge
    if (itPED == edges.end())throw MoleExcept(1040702,"Graph::delEdge","This graph doesn't contains the following edge : "+edge->toString());

    try
    {
        Vertex& ve1= edge->vertex1;
        ItEdge ited = find(ve1.links.begin(),ve1.links.end(),edge);
        if (ited == ve1.links.end()) throw MoleExcept(1030101,"Edge::~Edge","Given edge is not part of this vertex  \n"+ve1.toString()+"\n"+toString());
        ve1.links.erase(ited);
        Vertex& ve2= edge->vertex2;
        ited = find(ve2.links.begin(),ve2.links.end(),edge);
        if (ited == ve2.links.end()) throw MoleExcept(1030101,"Edge::~Edge","Given edge is not part of this vertex \n"+ve2.toString()+"\n"+toString());
        ve2.links.erase(ited);
    }
    catch (MoleExcept &e)
    {
        Moleaccess=Levels::FATAL;
        throw;
    }


// Deleting the edge from Graph :
    delete *itPED;
    edges.erase(itPED);


// Updating the maxNumEd - Unecessary if the edge num is not the maxNumEd
    if (NumEd+1 != maxNumEd) return;
    maxNumEd=1;
    bool zer=false;
    for(ItEdge it2 = edges.begin(); it2 != edges.end(); it2++)
      {
        if ((*it2)->getNum() >= maxNumEd) {maxNumEd= (*it2)->getNum();if (maxNumEd==1) zer=true;}
      }
    if (maxNumEd >1|| zer==true) maxNumEd++;
}


void Graph::delEdge(const EdgeList& EdList) 
{
    const size_t n=EdList.size();
    for (size_t i=0; i < n; i++)
    {
        const Edge* const edge=EdList.at(0);
        if (edge == (Edge*)NULL) throw MoleExcept(1040801,"Graph::delEdge","No edge given");
        // Check edge existence
        const ItEdge itPED=find(edges.begin(), edges.end(),edge); // Deletion if it finds in edges the corresponding edge
        if (itPED == edges.end())throw MoleExcept(1040802,"Graph::delEdge","This graph doesn't contains the following edge : "+edge->toString());

        try
        {
            Vertex& ve1= edge->vertex1;
            ItEdge ited = find(ve1.links.begin(),ve1.links.end(),edge);
            if (ited == ve1.links.end()) throw MoleExcept(1030101,"Edge::~Edge","Given edge is not part of this vertex  \n"+ve1.toString()+"\n"+toString());
            ve1.links.erase(ited);
            Vertex& ve2= edge->vertex2;
            ited = find(ve2.links.begin(),ve2.links.end(),edge);
            if (ited == ve2.links.end()) throw MoleExcept(1030101,"Edge::~Edge","Given edge is not part of this vertex \n"+ve2.toString()+"\n"+toString());
            ve2.links.erase(ited);
        }
        catch (MoleExcept &e)
        {
            Moleaccess=Levels::FATAL;
            throw;
        }

        // Deleting the edge from Graph :
        delete *itPED;
        edges.erase(itPED);

    }

// Updating the maxNumEd - Unecessary if the edge num is not the maxNumEd
    maxNumEd=1;
    bool zer=false;
    for(ItEdge it2 = edges.begin(); it2 != edges.end(); it2++)
      {
        if ((*it2)->getNum() >= maxNumEd) {maxNumEd= (*it2)->getNum();if (maxNumEd==1) zer=true;}
      }
    if (maxNumEd >1|| zer==true) maxNumEd++;
}

/*! \fn void   Graph::delEdge  (const Edge&  ed) 
  * \brief delete the given edge from this graph
  * \throw 1040301 - Given edge not part of this graph
  * \throw 1020101 when the given edge is not found in one of its two vertexs
  * \param ed : Edge to delete
  * Remove the given edge from the graph.
  * It will first search in the graph if this edge is associated to this graph, then delete the edge
  * The deletion of the edge automatically call delEdge() function of vertexs involved.
  * Also update the maxNumEd value
  */
void  Graph::delEdge(const       Edge &  edge) 
{
  const unsigned int NumEd= edge.getNum();

  // Check edge existence
  const ItEdge itPED=find(edges.begin(), edges.end(),&edge); // Deletion if it finds in edges the corresponding edge
  if (itPED == edges.end())throw MoleExcept(1040901,"Graph::delEdge","This graph doesn't contains the following edge : "+edge.toString());
  try
  {
      Vertex& ve1= edge.vertex1;
      ItEdge ited = find(ve1.links.begin(),ve1.links.end(),&edge);
      if (ited == ve1.links.end()) throw MoleExcept(1030101,"Edge::~Edge","Given edge is not part of this vertex  \n"+ve1.toString()+"\n"+toString());
      ve1.links.erase(ited);
      Vertex& ve2= edge.vertex2;
      ited = find(ve2.links.begin(),ve2.links.end(),&edge);
      if (ited == ve2.links.end()) throw MoleExcept(1030101,"Edge::~Edge","Given edge is not part of this vertex \n"+ve2.toString()+"\n"+toString());
      ve2.links.erase(ited);
  }
  catch (MoleExcept &e)
  {
      Moleaccess=Levels::FATAL;
      throw;
  }
  delete *itPED;
  edges.erase(itPED);

  // Updating the maxNumEd - Unecessary if the edge num is not the maxNumEd
  if (NumEd+1 != maxNumEd) return;
  maxNumEd=1;
  bool zer=false;
  for(ItEdge it2 = edges.begin(); it2 != edges.end(); it2++)
    {
      if ((*it2)->getNum() >= maxNumEd) {maxNumEd= (*it2)->getNum();if (maxNumEd==1) zer=true;}
    }
  if (maxNumEd >1|| zer==true) maxNumEd++;
}




/** \fn void Graph::renumEdge()
  \brief Renum all edges of the graph

  Scan all edges of the graph and renumerotate all. Call setNum() function for all edges
  */
void Graph::renumEdge()
{
  maxNumEd=0;
  for (ItEdge it= edges.begin(); it != edges.end(); it++){
      ((*it))->num=maxNumEd;
      maxNumEd++;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// MISCELLANEOUS ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/** \fn std::string Graph::toString() const
  \brief Return a text describing the graph

  Scan all vertexs and edges, call for each the toString() function and return a full string describing the graph
  */
std::string Graph::toString() const
{
  ostringstream oss;oss << "########## GRAPH "<< this<< "##########"<<endl;
  for (ItCVert it = vertexs.begin(); it != vertexs.end(); it++) oss << (*it)->toString()<<endl;
  for (ItCEdge it = edges.begin(); it != edges.end(); it++) oss <<"EDGE-"<< (*it)->toString()<<endl;
  return oss.str();
}




/*! \fn void   Graph::createCopyFrom(Graph *const Gr, bool with_parent)
    * Copy the given graph into this graph
    * \param Gr : Graph to copy from
    * \param with_parent : When cloning, also add related information for edge and vertex (see Clone() functions)
    * \throw MoleExcept Fatal error when unable to create a copy
    */
void   Graph::createCopyFrom(const Graph &Gr)
{
  try
  {
        const size_t VeSize=Gr.vertexs.size();
        size_t i=0;
        map<const Vertex*,const Vertex*> OldToNew;
        addVertexs(VeSize);
        for (; i< VeSize;++i)
        {
            vertexs.at(i)->clone(*Gr.vertexs.at(i));
            OldToNew.insert(pair<Vertex*,Vertex*>(Gr.vertexs.at(i),vertexs.at(i)));
        }

        for (i=0; i<Gr.edges.size();++i)
        {
            const Edge &GrEd=*Gr.edges.at(i);
             Vertex &ve= const_cast<Vertex&>(*OldToNew.at(&GrEd.getVertex1()));
             Vertex &ve2= const_cast<Vertex&>(*OldToNew.at(&GrEd.getVertex2()));
            addEdge(ve,ve2).clone(GrEd);
        }
  }
  catch (MoleExcept &e)
  {
    e.addTrace("Graph::createCopyFrom()");
    cerr << "Unable to create a copy of the given graph in Graph::createCopyFrom due to : "<<endl
         << e.getCode()<<"\t"<<e.getSource()<<"\t"<<e.getData()<<endl;
    e.setCode(999999);
    throw;
  }

}


double Graph::getDensity() const
{
    const double eds=(double)edges.size(), ves=(double)vertexs.size();
    return eds*2/(ves*ves-1);

}

const Edge& Graph::getEdge(const size_t &pos) const 
{
    if (pos >= edges.size())
        throw MoleExcept(1041001,
                         "Graph::getEdge",
                         "Given position is above the number of edges");
    return *edges.at(pos);

}





/**
 * @brief Create distance matrix
 * @param matrix:  distance matrix as mono dimensional array
 *
 * getDistMatrix generates the distance matrix of the calling graph.
 * To do so, it scan all vertexs looking for their linked vertexs and assign
 * to them a distance of 1. Then it goes into a loop to fetch linked vertexs
 * of theses linked vertexs and add a new level. It performs this until no vertexs
 * left.
 * @warning calls renumVertex function
 *
 */
void Graph::getDistMatrix(std::vector<unsigned short>& matrix)
{
    const  size_t sizeVertex=vertexs.size();


    if (sizeVertex < 2)     throw MoleExcept(8001,"GRIFPS::Create_Matrix", "Not enough pairs");

    // Prepare the given matrix by allocating the matrix size
    matrix.clear();
    matrix.reserve(sizeVertex*sizeVertex);
    // Filling it with dummy values (1000)
    matrix.assign(sizeVertex*sizeVertex,1000);

    // We nned to renumerotate vertexs :
    renumVertex();

    // vertexdone is a boolean array used to avoid redundance during graph scan
    vector<bool>vertexdone(sizeVertex,false);

    // vertextodo is the vector of vertex that we will currently look at
    // while vertextmp is the list of vertex that we will look at in the next loop
    vector<const Vertex*> vertextodo,vertextmp;
    // Distance value from reference vertex to the other vertex
    unsigned short level=1;


    // Scanning all vertexs :
    for (size_t iRef=0; iRef< sizeVertex;++iRef)
    {

         const Vertex &veRef = *vertexs.at(iRef);
         const unsigned int &posRef=veRef.getNum();
        // Vertex against itself has a distance of 0
         matrix.at(posRef*sizeVertex+posRef)=0;

         // initialize test :
         vertextodo.clear();
         level=0;
         for (size_t pos=0; pos < sizeVertex;++pos) vertexdone.at(pos)=false;
        vertexdone.at(posRef)=true;

        // Scanning all linked vertex from reference vertex :
        for (size_t iRefVe=0;iRefVe < veRef.numEdges();++iRefVe)
        {
            const Vertex &veLink=veRef.getVertex(iRefVe);
            // They have a distance of 1:
            matrix.at(posRef*sizeVertex+veLink.getNum())=1;
            matrix.at(veLink.getNum()*sizeVertex+posRef)=1;
            // Pushing them into vertextodo so we can look at their edges
            vertextodo.push_back(&veLink);
        }


        while (!vertextodo.empty())
        {
            //  vertex to look in the next loop
            vertextmp.clear();

            level++;

            for (size_t iCompVe=0; iCompVe < vertextodo.size();++iCompVe)
            {
                const Vertex &veComp = *vertextodo.at(iCompVe);
                const unsigned int &posComp= veComp.getNum();
                // If we already did them (it's a closer path), we skip it
                if (vertexdone.at(posComp)) continue;
                vertexdone.at(posComp)=true;

                // Setting up distance :
                matrix.at(posRef*sizeVertex+posComp)=level;
                matrix.at(posComp*sizeVertex+posRef)=level;

                // Fetching all edges of the comparison one:
                for (size_t iCompVe=0;iCompVe < veComp.numEdges();++iCompVe)
                {
                    const Vertex &veLink=veComp.getVertex(iCompVe);
                    if (vertexdone.at(veLink.getNum())) continue;
                    // And pushing them into vertextmp
                    vertextmp.push_back(&veLink);
                }
            }
            // Switching vertextmp to vertextodo
            vertextodo.clear();
            if (vertextmp.empty())break;
            vertextodo=vertextmp;

        }// END WHILE


    }//END iRef


}// END getDistMatrix



/**
 * @brief Delete all edges and vertex
 */
void Graph::clear()
{
    for(size_t i = 0; i < edges.size(); ++i)   delete edges[i];
    for(size_t i = 0; i < vertexs.size(); ++i) delete vertexs[i];
  edges.clear();
  vertexs.clear();
}

