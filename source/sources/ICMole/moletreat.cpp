#include <iostream>
#include <algorithm>
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/graph.h"
#include "headers/ICMole/bond.h"
#include "headers/ICMole/vertex.h"
#include "headers/ICMole/edge.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/box.h"
#include "headers/ICPars/molereader.h"
//#define ICHEM_DEBUG 1
#ifdef ICHEM_DEBUG
#include "headers/ICPars/molewriter.h"
#endif
using namespace std;
using namespace ICMole;



static const struct
{
    string Atm1;
    string Res1;
    string Atm2;
    string Res2;
    unsigned int BdType;

} listMod[76]=
{
{"C" ,  "ACT",  "O"     ,"ACT", 1},//ADDED JD-180214 - 4J4Q
{"C" ,  "ACT",  "OXT"   ,"ACT", 1},//ADDED JD-180214 - 4J4Q
{"C" ,  "ACT",  "CH3"   ,"ACT", 1},//ADDED JD-180214 - 4J4Q
{"C2",	"ACM",	"SG"	,"CYS",	1},
{"O1P",	"AMP",	"CZ"	,"TYR",	1},
{"P",	"AMP",	"NZ"	,"LYS",	1},
{"P",	"AMP",	"OH"	,"TYR",	1},
{"P",	"AMP",	"OH"	,"TYR",	1},
{"ND2",	"ASN",	"C1"	,"NAG",	1},//ADDED JD-180214 - 4J4Q
{"C1",	"ATO",	"SG"	,"CYS",	1},
{"C11",	"B59",	"NZ"	,"LYS",	0},
{"C1",	"B9D",	"OD1"	,"ASP",	0},
{"CAW",	"BF0",	"NZ"	,"LYS",	1},
{"C1",  "BMA",  "O4"    ,"NAG", 1},//ADDED JD-180214 - 4J4Q
{"O3",  "BMA",  "C1"    ,"MAN", 1},//ADDED JD-180214 - 4J4Q
{"C8",	"CES",	"OG"	,"SER",	1},
{"S1P",	"COA",	"SG"	,"CYS",	1},
{"SG",	"CYS",	"C6"	,"FAD",	1},
{"SG",	"CYS",	"SG2"	,"GSH",	1},
{"SG",	"CYS",	"C28"	,"QUM",	1},
{"SG",	"CYS",	"S1"	,"S2H",	1},
{"SG",	"CYS",	"C1"	,"PLM",	1},//ADDED JD-180214 - 4J4Q
{"C1",	"CYX",	"SG"	,"CYS",	1},
{"S1",	"D1T",	"SG"	,"CYS",	1},
{"C1",	"DIF",	"CE1"	,"HIS",	0},
{"C1",	"F4P",	"OG"	,"SER",	1},
{"C6",	"FAD",	"SG"	,"CYS",	1},
{"C8M",	"FAD",	"SG"	,"CYS",	1},
{"C8M",	"FAD",	"ND1"	,"HIS",	1},
{"C2",	"FCN",	"SG"	,"CYS",	1},
{"C1",	"FLH",	"NZ"	,"LYS",	1},
{"C1",	"FLH",	"CE1"	,"TYR",	1},
{"C1",	"G2F",	"OE1"	,"GLU",	1},
{"C1",	"G2F",	"OE2"	,"GLU",	1},
{"C1",	"G3H",	"SG"	,"CYS",	1},
{"C1",	"GER",	"SG"	,"CYS",	1},
{"SG2",	"GSH",	"SG"	,"CYS",	1},
{"SG2",	"GTT",	"SG"	,"CYS",	1},
{"SD",	"HCS",	"SG"	,"CYS",	1},
{"ND1",	"HIS",	"C8M"	,"FAD",	1},
{"NE2",	"HIS",	"C1"	,"FLH",	1},
{"CZ",	"LN5",	"SG"	,"CYS",	1},
{"CZ",	"LN7",	"SG"	,"CYS",	1},
{"NZ",	"LYS",	"C4A"	,"PLP",	2},
{"NZ",	"LYS",	"C15"	,"RET",	2},
{"C1",	"MAN",	"OG"	,"SER",	1},
{"S1",	"MTN",	"SG"	,"CYS",	1},
{"C1",  "NAG",  "O4"    ,"NAG", 1},//ADDED JD-180214 - 4J4Q
{"C6",	"NDU",	"SG"	,"CYS",	1},
{"NI",	"NI",	"NE2"	,"HIS",	0},
{"NI",	"NI",	"N"     ,"GLN",	0},
{"NI",	"NI",	"N"     ,"HIS",	0},
{"NI",	"NI",	"ND1"	,"HIS",	0},
{"NI",	"NI",	"N"     ,"THR",	0},
{"C8",	"NRO",	"SG"	,"CYS",	1},
{"C9",	"OCR",	"SG"	,"CYS",	1},
{"CE",	"ONL",	"SG"	,"CYS",	1},
{"C4A",	"PLP",	"NZ"	,"LYS",	2},
{"C10",	"PRF",	"SG"	,"CYS",	1},
{"N8",	"0T4",	"SG"	,"CYS",	1},
{"SD",	"2A2",	"SG"	,"CYS",	1},
{"C26",	"QUM",	"OE2"	,"GLU",	1},
{"C28",	"QUM",	"OD2"	,"ASP",	1},
{"C15",	"RET",	"NZ"	,"LYS",	2},
{"S1",	"S3H",	"SG"	,"CYS",	1},
{"C7",	"SA2",	"OG"	,"SER",	1},
{"OG",	"SER",	"C1"	,"MAN",	1},
{"C7",	"SFR",	"OG"	,"SER",	1},
{"OG1",	"THR",	"C1"	,"MAN",	1},
{"C6",	"UMP",	"SG"	,"CYS",	1},
{"C6'",	"UPG",	"SG"	,"CYS",	1},
{"C4A",	"X04",	"NZ"	,"LYS",	2},
{"C11",	"XLX",	"SG"	,"CYS",	1},
{"C1",	"ZPR",	"OG"	,"SER",	1},
{"C1",	"ZAH",	"OG"	,"SER",	1},
{"CS",	"XCN",	"SG"	,"XCN",	1}};



/**
 * @brief createMoleGraph
 * @param graph : Graph to fill in with molecular structure
 * @param onlyUsed : Takes only atom and bond with inUse set to true (default true)
 * @param wHydrogen : Set to true if hydrogen are included into the graph
 *
 *
 * Will create a graph from the molecule. Each vertex will be an atom and each
 * edge will be a bond.
 *
 */
void Molecule::createMoleGraph(Graph& graph,
                               const bool& onlyUsed,
                               const bool& wHydrogen) const
{

    try
    {

        Vertex* AtoV[maxNumAtom];
        for (size_t i=0; i<maxNumAtom;i++) AtoV[i]=(Vertex*)NULL;

        // Scanning each atoms :
        for (size_t i=0; i<Atoms.size();i++)
        {
            Atom* atm = Atoms.at(i);

            // Continue when unwanted
            if ((onlyUsed   && !atm->isUsed()) ||
                    (!wHydrogen && atm->isHydrogen()))  continue;

            // Create the vertex
            Vertex&ve = graph.addVertex();
            ve.setAtom(atm);
            AtoV[atm->getNum()]=&ve;
        }

        for (ItCBond it = Bonds.begin(); it != Bonds.end(); it++)
        {
            Bond* bd=*it;
            if (AtoV[bd->getAtom1().getNum()] == (Vertex*)NULL
                    ||AtoV[bd->getAtom2().getNum()] == (Vertex*)NULL)continue;
            Vertex &ve1 = *(AtoV[bd->getAtom1().getNum()]),
                    &ve2 = *(AtoV[bd->getAtom2().getNum()]);
            // cout << ve1.toString()<< ve2.toString()<<endl;
            graph.addEdge(ve1,ve2
                          ).setBond(bd);




        }
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Molecule::createMoleGraph");
        throw;

    }

}

#ifdef ICHEM_DEBUG
void saveM( Graph& graph, const unsigned int N)
{
    Molecule mole;
    mole.setMoleType(MoleType::LIGAND);
    map<Vertex*,Atom*> OtoN;
    for (ItCVert it = graph.firstVertex(); it != graph.lastVertex();++it)
    {
        Vertex* ve=*it;
        Atom &atm=mole.addAtom(ve->getAtom()->getAtomicName(),
                               ve->getAtom()->fixpos,
                               ve->getAtom()->getName(),
                               ve->getAtom()->getMOL2Type(),
                               &mole.addResidu(ve->getAtom()->getResidu()->getChainName(),ve->getAtom()->getResidu()->getFNum(),ve->getAtom()->getResidu()->getName()));
        OtoN.insert(pair<Vertex*,Atom*>(ve,&atm));
    }
    for (ItCEdge it = graph.firstEdge(); it != graph.lastEdge();++it)
    {
        Edge* ed=*it;
        mole.addBond(OtoN.at(&ed->getVertex1()),OtoN.at(&ed->getVertex2()),BondType::SINGLE);

    }
    MoleWriter mw;
    ostringstream oss;
    oss<< "TEST_"<<N<<".mol2";
    mw.newFile(oss.str());
    mw.writeMOL2(&mole);
}
#endif












/**
 * @brief Molecule::ringPerception
 *
 *
 * Detect all cycles within the given molecule. This is made by generating a
 * molecular graph and reducing it to fetch cycles.
 * Cycle object will be created for each cycle and stored into molecule.
 * Moreover, atomic properties will be set to tell whether the atom is within
 * a scaffold, a substituent or a linker.
 *
 */
void Molecule::ringPerception() throw(MoleExcept)
{
    if (typemol == MoleType::PROTEIN)
    {
        const string AA_AROM =  " PHE TYR HIS TRP dA dG rA rG dC dT rC rU ";
        string list;
        AtomList alist;

        for (ItCRes itR = Residues.begin(); itR != Residues.end(); itR++)
            {
            Residu *res = *itR;
//            cout << res->getIdentifier() << endl;
                if (AA_AROM.find(res->getName().substr(0,3)) == string::npos){ continue;}
                if (res->getName().find("PHE") != string::npos) list=" CD1 CE1 CZ CE2 CD2 CG ";
                else if (res->getName().find("TYR") != string::npos) list=" CD1 CE1 CZ CE2 CD2 CG ";
                else if (res->getName().find("HIS") != string::npos) list=" CG CD2 ND1 CE1 NE2 ";
                else if (res->getName().find("TRP") != string::npos) list=" CD2 CE2 CZ2 CE3 CH2 CZ3 ";
                else if (res->getName().find("dA dG rA rG") != string::npos) list=" N1 C2 N3 C4 N9 C8 N7 C5 C6 ";
                else if (res->getName().find("dC dT rC rU") != string::npos) list=" N1 C2 N3 C4 C5 C6 ";
                alist.clear();

                for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom();itA++)
                    {
                    Atom *atm = *itA;
                        if (list.find(" "+atm->getName()+" ") == string::npos) continue;
                        alist.push_back(atm);

                    }

                if (alist.size() < 3) continue;
                addCycle(alist);
            }
        return;
    }



////////////////////////////////////////////////////////////
//////////Step 1 : Create the graph without hydrogen :
////////////////////////////////////////////////////////////
    Graph graph(Atoms.size(),Bonds.size());
    createMoleGraph(graph,true,false);


    // STEP 2 : Create an atom to vertex mapping :
    vector<Atom*> VtoA;
    for (size_t i=0;i<graph.numVertex();i++)
    {
        const Vertex &ve = graph.getVertex(i);
        VtoA.push_back(ve.getAtom());
    }


////////////////////////////////////////////////////////////
////////// Step 2 : Delete all vertex that have only one edge :
////////////////////////////////////////////////////////////
    bool changes=true;
// List of vertex that will be deleted
    VertexList toDel;

    // While we find a vertex to delete, we continue
    while(changes)
    {
        changes=false;
        toDel.clear();
        // Scanning each vertex of the graph:
        for (ItVert it  = graph.firstVertex();
             it != graph.lastVertex();
             ++it)
        {
            Vertex &ve = **it;
            // To check if they have 1 edge
            if (ve.numEdges() > 1)continue;

#ifdef ICHEM_DEBUG
            cout << "Delete : "<< ve.getAtom()->getIdentifier()<<endl;
#endif
            // If so, it is a substituent and we delete the vertex:
            ve.getAtom()->props.setSubstituent(true);
            toDel.push_back(&ve);

        }
        // When no vertex can be found we stop:
        if (toDel.empty())break;
        changes=true;
#ifdef ICHEM_DEBUG
        cout << "Graph size before : "<< graph.numVertex()<<endl;
#endif
        graph.delVertexs(toDel);
#ifdef ICHEM_DEBUG
        cout << "Graph size after : "<< graph.numVertex()<<endl;
#endif
    }// END WHILE CHANGE


////////////////////////////////////////////////////////////
////////// Step 3 : Assign to each edge its involved vertex:
////////////////////////////////////////////////////////////
    ostringstream oss;
    for (ItEdge it = graph.firstEdge();it != graph.lastEdge();it++)
    {
        oss.str("");
        oss <<  (*it)->getVertex1().getNum()<<"-"
            <<  (*it)->getVertex2().getNum();
        (*it)->setLabel(oss.str());
    }


////////////////////////////////////////////////////////////
////////// Step 4 : Reduce graph
////////////////////////////////////////////////////////////
/// Principle:
///    When we have a vertex with two edges (e.g. VE2) : VE1-VE2-VE3-VE4-VE5
///    With VE5 linked to VE1
///    We can remove VE2 and create a new edge ED1 between VE1 and VE3
///    ED1 has a label that will be set to VE2
///    So we have smt like that: VE1--(ED1)---VE3-VE4-VE5 with ED1=>VE2
///    When we continue, we remove VE3 and create ED2 between VE1 and VE4
///    So we have : VE1---(ED2)---VE4-VE5 with ED2=>VE2/VE3
///    Now we have a triangle since VE1 is linked to VE5
///    When a vertex (VE1) is linked to two other vertexs (VE4 and VE5)
///    and when theses two other vertexs are linked to each other
///    Then we have a cycle
///    To get all atoms within the cycle, we fetch label of the 3 edges
///    ED2, VE1-VE5 and VE4-VE5
///    So we get from labels : VE2/VE3
///    And the 3 vertexs involved: VE1/VE4/VE5
///    Therefore we have our cycle and all atoms involved in the cycle
///    We the remove VE1 leading to VE4-VE5
///    Since they have 1 edge, they are both deleted and the process is over.

    bool alldone;
    vector<unsigned int> ListAt;
    // 3 steps is more than enough to finish the process:
    for (unsigned int IVar1=0; IVar1 < 4; IVar1++)
    {
        alldone=false;
        while (!alldone)// Continue until no change
        {
            alldone = true;
            // Scanning each vertex:
            for (size_t ive=0;ive<graph.numVertex();ive++)
            {
                Vertex& ve= const_cast<Vertex&>(graph.getVertex(ive));
#ifdef ICHEM_DEBUG
                cout << "LOOKING AT : ("<<ve.numEdges()<< " edges) "<<ve.toString()<<endl;
#endif
                // Delete the vertex is one edge or none
                if (ve.numEdges() <= 1) {
#ifdef ICHEM_DEBUG
                    cout << "DELETING - " << ve.numEdges() <<" edges"<<endl;
#endif
                    graph.delVertex(&ve,false);alldone=false;
                }
                // When two edges:
                else if (ve.numEdges() == 2)
                {
                      Vertex&ve1= ve.getVertex(0);
                      Vertex&ve2= ve.getVertex(1);
                    const Edge *ed1 = ve.getEdge(0);
                    const Edge *ed2 = ve.getEdge(1);
                    // Check that the two other vertex are linked to each other
                    if (ve1.hasEdgeWith(ve2))
                    {

                        // Therefore Cycle:
                        // We fetch labels:
                        oss.str("");
                        oss << ed1->getLabel() <<"-"
                            << ed2->getLabel()<<"-"
                            << ve1.getEdgeWith(ve2)->getLabel();
                        ListAt.clear();
                        ICMole::tokenize(oss.str(),ListAt,"-");
                        // And convert it to atoms:
                        AtomList ATML;
                        for (vector<unsigned int>::iterator it = ListAt.begin();
                                                            it != ListAt.end();
                                                            it++)
                        {
                            ATML.push_back(VtoA.at((*it)));
                        }
                        // And add cycle:
                        addCycle(ATML);

                    }
                    // Otherwise, we create a new edge
                    // between the two other vertex
                    // And remove the former one
                    else
                    {
                        Edge &new_ed = graph.addEdge(ve1,ve2);// Creating a new edge

                        oss.str("");
                        oss<< ed1->getLabel()<<"-"
                           << ed2->getLabel();
                        new_ed.setLabel(oss.str());// Setting label
#ifdef ICHEM_DEBUG
                        cout << "DELETING : "<< ve.toString(true)<< " => New edge : " <<new_ed.toString()<<endl ;
#endif

                    }
                    graph.delVertex(&ve,false);
                    alldone=false;


                }

            }
        }
        // When no reduction can be made from 0,1,2 edges
        // It can be because we have only vertex with 3 edges
        alldone=false;
        while (!alldone)
        {
            alldone=true;
            for (size_t ive=0;ive<graph.numVertex();ive++)

            {
                Vertex& ve= const_cast<Vertex&>(graph.getVertex(ive));
                if (ve.numEdges()!=3)continue;

                // Scanning all edges to find
                // Two other vertexs linked to each other:
                for (size_t I1=0; I1 <3; I1++)
                {
                    Vertex &ve1 = ve.getVertex(I1);
                    const Edge* ed1 = ve.getEdge(I1);
                    for (size_t I2=I1+1; I2 <3; I2++)
                    {
                        Vertex &ve2 = ve.getVertex(I2);
                        const Edge* ed2 = ve.getEdge(I2);
                        if (ve1.hasEdgeWith(ve2))
                        {

                            // CYCLE
                            oss.str("");
                            oss << ed1->getLabel() <<"-"
                                << ed2->getLabel() <<"-"
                                << ve1.getEdgeWith(ve2)->getLabel();
                            ListAt.clear();
                            ICMole::tokenize(oss.str(),ListAt,"-");
                            AtomList ATML;
                            for (vector<unsigned int>::iterator
                                 it = ListAt.begin();it != ListAt.end(); it++)
                            {
                                ATML.push_back(VtoA.at((*it)));
                            }
                            addCycle(ATML);
                        }
                        else
                        {
                            Edge &new_ed = graph.addEdge(ve1,ve2);// Creating a new edge

                            oss.str(""); oss<< ed1->getLabel()<<"-"
                                            <<ed2->getLabel();
                            new_ed.setLabel(oss.str());// Setting label
#ifdef ICHEM_DEBUG
                            cout <<  " => New edge : " <<new_ed.toString()<<endl ;
#endif
                        }
                    }// END size_t I2
                }// END size_t I1

                alldone =false;

                graph.delVertex(&ve,false);


            }
        }
    }
#ifdef ICHEM_DEBUG
    cout <<"-----------------"<<endl<<"-----------------"<<endl<<"-----------------"<<endl;
    cout << graph.toString();
#endif


    if (graph.numEdges() >= 1)    cerr << "WARNING : Unable to find all cycles within the molecule "<< name<<endl;



  for (ItCAtom itA = Atoms.begin(); itA != Atoms.end(); itA++)
  {
      Atom &atmA = **itA;
      if (atmA.isHydrogen())continue;
      if (!atmA.props.isScaffold() && !atmA.props.isSubstituent())
          atmA.props.setLinker(true);
  }
  for (ItCAtom itA = Atoms.begin(); itA != Atoms.end(); itA++)
  {
      Atom &atmA = **itA;
      if (!atmA.isHydrogen() || atmA.props.isDummy())continue;
      const Atom &atmB = atmA.getAtomLinked(0);
      if (atmB.props.isScaffold()) atmA.props.setScaffold(true);
      else if (atmB.props.isLinker())atmA.props.setLinker(true);
      else if (atmB.props.isSubstituent())atmA.props.setSubstituent(true);

  }


}












void Molecule::addCycle(const AtomList &ListAt)
{

    size_t NbrDouble=0, NbrArom=0;
    const size_t NbrAtm = ListAt.size();
    Cycle *cycle= new Cycle(this);
    cycle->addAtoms(ListAt);
    BondList listBd;

#ifdef ICHEM_DEBUG
    cout << "NEW CYCLE : "<<endl;
#endif
    for (ItCAtom itA = ListAt.begin();
         itA != ListAt.end();
         itA++)
    {
        Atom &atmA = **itA;
#ifdef ICHEM_DEBUG
        cout << atmA.getIdentifier()<<endl;
#endif
        atmA.props.setRing(true);
        atmA.props.setScaffold(true);
        atmA.inCycle=true;
        ItCAtom itB=itA;
        for (++itB; itB != ListAt.end(); itB++)
        {
            Atom &atmB = **itB;
            if (!atmA.hasBondWith(atmB))continue;
            Bond *bd=const_cast<Bond*>(atmA.getBondWith(atmB));
            bd->props.setRing(true);
            if (bd->getBondType() == BondType::DOUBLE)   NbrDouble++;
            else if (bd->getBondType() == BondType::AROMATIC) NbrArom++;
            listBd.push_back(bd);
        }
    }
#ifdef ICHEM_DEBUG
    cout << NbrAtm<< " " << NbrDouble << " " << NbrArom<<endl;
#endif
    if ((NbrAtm == 6 && (NbrDouble==3|| NbrArom ==6))// CHECK
            || (NbrAtm == 6 && (NbrDouble==4 || NbrArom == 8))
            || (NbrAtm == 4 && (NbrDouble==4 || NbrArom == 8))
            || (NbrAtm == 5 && (NbrDouble==2 && NbrArom == 2))
            || (NbrAtm == 5 &&  NbrDouble==2 && NbrArom == 0)// CHECK
            || (NbrAtm == 5 &&  NbrDouble==0 && NbrArom ==10))
    {
        cycle->setAromatic(true);
        aromaticRes.addAtom(cycle->getCenter());
        for (ItBond it = listBd.begin(); it != listBd.end(); it++)
            (*it)->props.setAromatic(true);
#ifdef ICHEM_DEBUG
        cout << "AROMATIC"<<endl
             << "END CYCLE"<<endl;
#endif
    }
    else
    {
        cycle->getCenter().setName("DuCy");
        cyclicRes.addAtom(cycle->getCenter());
#ifdef ICHEM_DEBUG
        cout << "ALYPHATIC"<<endl
             << "END CYCLE"<<endl;
#endif
    }
    cycle->getFixpos();
    Cycles.push_back(cycle);
}



void Molecule::checkMOL2() throw(MoleExcept)
{

  for (ItCAtom itA = Atoms.begin(); itA != Atoms.end(); ++itA)
    {try
      {
          (*itA)->checkMOL2type();
      }catch (MoleExcept &e)
      {
        e.addTrace("Molecule::checkMOL2");
        switch (ICMole::Moleaccess)
          {
          //case Levels::QUIET:break;
          case Levels::NOTICE:cerr<< e.getData()<<endl;
            cerr << (*itA)->toString()<<endl;
            break;
          case Levels::STRICT:
          case Levels::FATAL:throw;
          }
      }

    }

}


void Molecule::clear()
{
    for(size_t i = 0; i < Cycles.size(); ++i) delete Cycles[i];
    for(size_t i = 0; i < Residues.size(); ++i) delete Residues[i];
    for(size_t i = 0; i < Chains.size(); ++i) delete Chains[i];
    for(size_t i = 0; i < Bonds.size(); ++i) delete Bonds[i];
    for(size_t i = 0; i < Atoms.size(); ++i) delete Atoms[i];
    Cycles.clear();
    Residues.clear();;
    Chains.clear();;
    Bonds.clear();;
    Atoms.clear();
    clearSets();
    maxNumAtom=0;
    maxNumBond=0;
    maxNumResidu=0;
    bary_check=false;

}


void Molecule::selChains(const std::vector<std::string>& allowedchains,
                         const bool &delete_all,
                         const std::string& verbose
                         ) throw(MoleExcept)
{
    const bool for_scPDB=!verbose.empty();
    if (for_scPDB) cout << verbose<<"|START\t"<<endl;

    AtomList listAtmToDel;

    for (ItCChain itChain = Chains.begin(); itChain != Chains.end(); ++itChain)
    {
        Chain &chain = **itChain;
        if (find(allowedchains.begin(), allowedchains.end(),chain.getName())
           != allowedchains.end())continue;
        if (for_scPDB) cout << verbose<<"|REMOVING\t"<< chain.getName()<<"\t";
        if (!delete_all)
        {
            chain.setUse(false,true);

        }
        else
        for (ItCRes itRes = chain.firstC(); itRes != chain.lastC(); ++itRes)
        {
           Residu &residu=**itRes;
           for (ItCAtom itAtom = residu.firstAtom();
                itAtom != residu.lastAtom();
                ++itAtom)
               listAtmToDel.push_back(*itAtom);
        }
        if (for_scPDB) cout << "END"<<endl;

    }
    if (delete_all && !listAtmToDel.empty()) delAtoms(listAtmToDel);



if (for_scPDB) cout << verbose<<"|END"<<endl;
}


/**
 * @brief Molecule::createConnect
 * @param cleanAll : Delete all bonds before running it
 * @param verbose : For log purpose, you can set up a string that will
 * be the header of all log lines
 *
 * Recreate connectivity by looking at VdW radius of each atoms
 *
 */
void Molecule::createConnect(ResiduList &failedRes, AtomList &failedAtom, const bool& cleanAll,
                             const string& verbose)  throw(MoleExcept)
{
    const bool wVerbose= !verbose.empty();
    bool useGrid=false;
    Grid *grid=(Grid*)NULL;
    if (complex != (Complex*)NULL && complex->uptoGrid)
    {
        grid=&complex->getGrid(0);
        useGrid=true;

    }
// Deletion of all existing bond.
  if (cleanAll) cleanBond();

  if (wVerbose){ cout << verbose<<"START\t"<<endl;}


  double dist=0, VdW=0;
  AtomList toComp;
  bool  atmBox=false;
  AtomList atomlist;

  for (size_t posR = 0; posR < Atoms.size(); ++posR)
  {
      Atom &atomR = *Atoms.at(posR);
      toComp.clear();
      atmBox=true;

      if (useGrid)
      {

          Box *box=atomR.getBox(grid);
          if (box == (Box*)NULL)
          {
              atmBox=false;
          }
          else
          {
              grid->getAdjacentAtoms(atomlist,atomR,5,true);

                for (ItCAtom itA = atomlist.begin(); itA != atomlist.end();++itA)
                {
                    Atom &atomC = **itA;

                    if (atomR.getNum() >=atomC.getNum()) continue;
                    toComp.push_back(&atomC);
                }
          }

      }
      if (!useGrid || !atmBox)
      {
          for (size_t posC = posR+1; posC < Atoms.size(); ++posC)
          {
              Atom &atomC = *Atoms.at(posC);
              if (atomC.fixpos.calcDist(atomR.fixpos,4.1) > 4) continue;
              toComp.push_back(&atomC);
          }
      }

      sortAndUnique(toComp);

size_t Nmade=0;
      const size_t nToComp = toComp.size();
      for (size_t posC=0; posC < nToComp; posC++)
      {
          Atom &atomC = *toComp.at(posC);

          // RULE 1 : 1 water molecule can only be link with itself
          if ((atomR.getResidu()->getResType()==ResType::WATER
            ||atomC.getResidu()->getResType()==ResType::WATER)
            &&atomR.getResidu() != atomC.getResidu()) continue;

          // RULE 2 : A metallic atom can only be linked with atoms of the same residu
          if ((atomR.isMetallic() && atomR.getResidu() != atomC.getResidu())
            ||(atomC.isMetallic() && atomR.getResidu() != atomC.getResidu())) continue;
          // 1 Hydrogen cannot be linked to an atom coming from another residu
          if ((atomR.isHydrogen() || atomC.isHydrogen())
             && atomR.getResidu() != atomC.getResidu()) continue;

          // A ionic molecule cannot be linked to another residu
          if ((atomR.getResidu()->getResType() == ResType::ION
            && atomC.getResidu() != atomR.getResidu()) ||
              (atomC.getResidu()->getResType() == ResType::ION
           && atomR.getResidu() != atomC.getResidu())) continue;

          // FeS cluster managment
          if (atomR.getAtomicName() == "Fe"     && atomC.getAtomicName()=="S") VdW = 2.5;
          else if (atomR.getAtomicName() == "S" && atomC.getAtomicName()=="Fe") VdW = 2.5;
          else VdW= (atomR.getVdWRadius()+atomC.getVdWRadius())/2+0.3;

          dist = atomR.fixpos.calcDist(atomC.fixpos,VdW+0.1) ;
          if (dist < VdW/2 &&
             (!atomR.isHydrogen() || !atomC.isHydrogen()) &&
             (atomR.getName().substr(0,1)!="H" && atomC.getName().substr(0,1)!="H"))
          {
              failedRes.push_back(atomR.getResidu());
              failedRes.push_back(atomC.getResidu());
              failedAtom.push_back(&atomR);
              failedAtom.push_back(&atomC);
              cerr << verbose<<"ERROR\tATOMS TOO CLOSE : "
                   <<atomR.getIdentifier()<<"\t"
                  <<atomC.getIdentifier()<<"\t"<<dist<<endl;
          }
          if (dist >= VdW) continue;
          Nmade++;
              try{
                if (atomR.isHydrogen())
                  {
                    if (atomR.getNumBond()==0)
                    {
                        addBond(atomR,atomC,BondType::SINGLE);
                    }
                    else
                      {
                        Bond *bd = atomR.getBond(0);
                        if (dist < bd->getLength())
                          {
                            delBond(bd);
                            addBond(atomR,atomC,BondType::SINGLE);
                          }

                      }
                  }
                else if (atomC.isHydrogen())
                  {
                    if (atomC.getNumBond()==0) addBond(atomR,atomC,BondType::SINGLE);
                    else
                      {
                        Bond *bd = atomC.getBond(0);
                        if (dist < bd->getLength())
                          {
                            delBond(bd);
                            addBond(atomR,atomC,BondType::SINGLE);
                          }

                      }
                  }
                else addBond(atomR,atomC,BondType::DUMMY);}

              catch (MoleExcept &e)
              {
                if (e.getCode() != 1060104){
                    cerr << verbose<<"|ERRORS\t"<<e.getCode()<<"\t"<<e.getData()<<endl;}
              }



      }
      //cout << Nmade<<endl;


  }


  // Correcting some knowledge based issues
  std::string RName1;
  const size_t nResidu = Residues.size();

  for (size_t nRes = 0; nRes < nResidu; ++nRes)
  {
      Residu &residu = *Residues.at(nRes);
      RName1 = residu.getName().substr(0,3);
      for (size_t Npos = 0; Npos < 68; ++Npos)
      {
          if (listMod[Npos].Res1 != RName1)continue;
          for (ItCAtom itAtom = residu.firstCAtom();
                       itAtom != residu.lastCAtom();
               ++itAtom)
          {
              Atom &atom1 = **itAtom;
              if (atom1.getName()==listMod[Npos].Atm1)continue;
              for (size_t NBond = 0; NBond < atom1.getNumBond();++NBond)
              {
                  Bond &bond = *atom1.getBond(NBond);
                  Atom &atom2= bond.getOtherAtom(atom1);
                  if (atom2.getName() != listMod[Npos].Atm2
                  ||  atom2.getResiduName() != listMod[Npos].Res2) continue;
                  switch (listMod[Npos].BdType)
                  {
                  case 0: delBond(&bond);  break;
                  case 1: bond.setBondType(BondType::SINGLE);break;
                  case 2 :bond.setBondType(BondType::DOUBLE);break;
                  }
                  break;
              }
          }

      }
  }




  if (wVerbose){ cout << verbose<<"END\t"<<endl;}

}

void Molecule::setResiduTypes(const string &verbose) throw(MoleExcept)
{
    const bool wVerbose= !verbose.empty();
    if (wVerbose)cout << verbose<<"START"<<endl;

    const size_t nResidu = Residues.size();
    for (size_t nRes = 0; nRes < nResidu;++nRes)
    {
        Residu& residu = *Residues.at(nRes);
        try
        {
        residu.applyResiduType(verbose);
        }catch (MoleExcept &e)
        {

                e.addTrace("Molecule::setResiduTypes");
                throw;



        }
    }
    if (wVerbose)  cout << verbose<<"END"<<endl;
}






void Molecule::matchTemplate(ResiduList& failedRes, AtomList& failedAtom, const string& verbose,const bool& applyMOL2type) throw(MoleExcept)
{

    const bool wVerbose= !verbose.empty();
    if (wVerbose)cout << verbose<<"START"<<endl;

    ////////// PARAMS //////////
// DATA RELATED TO TEMPLATE MANAGEMENT
    const std::string ICPath = get_IChem_LIB_Path();
    // Load HETLIST file to get all residues:
    Residu::loadHETClass();

    map<string,HetData>::iterator itHET;


// DATA RELATED TO MANAGEMENT OF 1 TEMPLATE:

    MoleReader              mread;
    struct templateData
    {
     Molecule mole;
     std::map<const Atom*,unsigned short> order;
     std::vector<unsigned short> distance;
     unsigned short nHeavy;

    } ;

// DATA RELATED TO RESIDU :
    const size_t nResidu = Residues.size();
    std::map<const Atom*,unsigned short> residuOrder;
    std::vector<unsigned short> residuDistance;


// DATA RELATED TO GRAPH MATCHING:

    CliqueData cliqueData;
    std::string AtmDone;


// STEP 1 - Getting list of all residu name :
    vector<string> listResidu;
    for (size_t nRes = 0; nRes < nResidu;++nRes)
    {
        const std::string &residuName=Residues.at(nRes)->getName();
        if (residuName == "HOH" || Residues.at(nRes)->getResType()==ResType::METAL)continue;
        listResidu.push_back(residuName);
    }
    // Sorting the list and make it unique
    // So we can have the exact number of templates to load:
    sort(listResidu.begin(),listResidu.end());
    vector<string>::iterator it=unique(listResidu.begin(),listResidu.end());
    listResidu.resize( std::distance(listResidu.begin(),it) );

    // templateCount is therefore the count of unique residues in the protein:
    const size_t templateCount=listResidu.size();
    // Creating an array of templateData structure
    templateData *templatesTable = new templateData[templateCount] ;
try{
// STEP 1 - Loading all templates
    for (size_t posTpl = 0; posTpl < templateCount;posTpl++)
    {
        try
        {
        const std::string& residuName = listResidu.at(posTpl);

        //cout << residuName<<" " << posTpl<<endl;
        templateData& templatedata= templatesTable[posTpl];
        itHET = Residu::HETClass.find(residuName);
        if (itHET == Residu::HETClass.end())
        {
            throw MoleExcept(4001011,
                             "MatchTemplate",
                             "TEMPLA\t"+residuName+" Cannot find HET code information" );
        }
        // If so, we can load the template:
        HetData &HD = (*itHET).second;
        if (!HD.hasTemplate) {
            throw MoleExcept(4001010,
                             "MatchTemplate",
                             "ERROR\tCannot find template for HET "+residuName);
          }
        // Loading the template:
        mread.loadNewFile(ICPath+"/datas/templates/"+residuName+".mol2");
        Molecule &tplmole=templatedata.mole;
        mread.loadNextMolecule(tplmole,MoleType::LIGAND);
        // Creating the molecular distance matrix and order:
        tplmole.checkMOL2();
        tplmole.renumAtom();
        if (tplmole.numResidus()==0)
        {
            Residu &newres=tplmole.addResidu("X",1,residuName);
            for (ItCAtom itA = tplmole.firstAtom();itA != tplmole.lastAtom();++itA) tplmole.setResToAtom(*itA,&newres);
        }
        templatedata.mole.getResidu(0).getDistMatrix(templatedata.distance,templatedata.order);
        templatedata.nHeavy=templatedata.order.size();
        }
        catch (MoleExcept &e)
        {
            cerr << e.getCode()<<endl
                 << e.getData()<<endl
                 << e.getSource()<<endl
                 << e.getTrace()<<endl;

        }
    }


// STEP 3 - Apply templates
    const size_t nResidu = Residues.size();

    for (size_t nRes = 0; nRes < nResidu;++nRes)
    {

        Residu& residu = *Residues.at(nRes);
        const std::string& residuName = residu.getName();

        const size_t nAtoms= residu.numAtom();
        try{
       // cout << residu.getIdentifier()<<endl;
// SIMPLE CASE => WATER :
        if (residu.getResType()==ResType::WATER)
        {
          //  cout << residu.getIdentifier()<< " WATER"<<endl;
            for (size_t nAtom =0; nAtom < nAtoms;++nAtom)
            {
                Atom& atom = residu.getAtom(nAtom);
                if (atom.isOxygen())atom.setMOL2Type("O.3");
                else atom.setMOL2Type("H");
            }
            continue;
        }
// SIMPLE CASE => METAL:
        else if (residu.getResType()==ResType::METAL)
        {
            if (nAtoms != 1) throw MoleExcept(10000,"Molecule::MatchTemplate",
                                             "Metal residu must contains 1 atom "+residu.getIdentifier());
            Atom& atom=residu.getAtom(0);
            if (atom.getName()=="CO")atom.setMOL2Type("Co.oh");
            else atom.setMOL2Type(atom.getAtomicName());
            continue;
        }
// Getting position of the residu in the template list:
        const size_t templatePos=std::distance(listResidu.begin(),std::find(listResidu.begin(),listResidu.end(),residuName));
// Getting template data:
        templateData& templatedata= templatesTable[templatePos];

// Getting residu distance matrix:

        residu.getDistMatrix(residuDistance,residuOrder);

AtmDone="";vector<const Atom*> atN;
        cliqueData.listPairs.clear();
        cliqueData.graph.clear();
        cliqueData.cliques.clear();
        for(std::map<const Atom*,unsigned short>::const_iterator
            itResAtom=residuOrder.begin();
            itResAtom!=residuOrder.end();
            ++itResAtom)
        {
            const Atom& atomRes=*(*itResAtom).first;
            atN.push_back(&atomRes);
            AtmDone+="0";
#ifdef ICHEM_DEBUG
            cout << "RESIDU:"<<atomRes.getName()<<" \t "<<(*itResAtom).second<<endl;
#endif
            for(std::map<const Atom*,unsigned short>::const_iterator
                itTplAtom=templatedata.order.begin();
                itTplAtom!=templatedata.order.end();
                ++itTplAtom)
            {
                const Atom& atomTpl=*(*itTplAtom).first;
                if (atomRes.getAtomicName()!=atomTpl.getAtomicName())continue;


                Pair paire(&atomRes,
                         &atomTpl,
                         &cliqueData.graph.addVertex());

                cliqueData.listPairs.push_back(paire);

            }
        }
#ifdef ICHEM_DEBUG
        for(std::map<const Atom*,unsigned short>::const_iterator
            itTplAtom=templatedata.order.begin();
            itTplAtom!=templatedata.order.end();
            ++itTplAtom)
        {
            const Atom& atomTpl=*(*itTplAtom).first;
            cout << "TEMPLATE:"<<atomTpl.getName()<<" \t "<<(*itTplAtom).second<<endl;
        }

#endif

        for (vector<Pair>::iterator itR = cliqueData.listPairs.begin(); itR != cliqueData.listPairs.end(); itR++)
        {
              const Pair& PeR = *itR;
              const unsigned short& rrefpos=residuOrder.at(PeR.ref);
              const unsigned short& rcomppos=templatedata.order.at(PeR.comp);
            for (vector<Pair>::iterator itC = itR+1; itC != cliqueData.listPairs.end(); itC++)
            {
                 Pair& PeC = *itC;

                  if (PeR.ref==PeC.ref || PeR.comp==PeC.comp)continue;
                 const unsigned short& crefpos=residuOrder.at(PeC.ref);
                 const unsigned short& ccomppos=templatedata.order.at(PeC.comp);
                 const unsigned short distRef = residuDistance.at(rrefpos*residuOrder.size()+crefpos);
                 const unsigned short distComp = templatedata.distance.at(rcomppos*templatedata.order.size()+ccomppos);
                 if (distRef!= distComp)continue;

                // cout << PeR.ref->getName()<<"<->"<<PeR.comp->getName()<<"\t"<<PeC.ref->getName()<<"<->"<<PeC.comp->getName()<<" " << distRef << " " << distComp<<endl;

                 cliqueData.graph.addEdge(*PeR.vertex,*PeC.vertex);


            }}

        if (cliqueData.listPairs.size() < 3) continue;
        GraphMatch GrMMatch(cliqueData);
        GrMMatch.setMinSizeClique(3);
        GrMMatch.createMatrix();

        GrMMatch.runSearchAllMaxCliques(false);

         if (cliqueData.cliques.size()==0) continue;

         unsigned int maxsize=0;
         vector<int > bestClic;
         for (vector<vector<int> >::const_iterator itCl = cliqueData.cliques.begin(); itCl != cliqueData.cliques.end(); itCl++)
         {
             if ((*itCl).size() <= maxsize) continue;
                 bestClic = *itCl;
                 maxsize = (unsigned int)(*itCl).size();
         }

         unsigned short nAtChanged=0;
         for (vector<int >::const_iterator itMat = bestClic.begin(); itMat!= bestClic.end(); itMat++)
         {
             Pair &Pe = cliqueData.listPairs.at(*itMat);
             nAtChanged++;
            // if (wVerbose)cout << verbose<<"MODIFS\t"<<residu.getIdentifier()<<"\tATOM\t"<< Pe.ref->getName()<<"\t"<<Pe.comp->getName()<<"\t"<<Pe.comp->getMOL2Type()<<endl;
           #ifdef ICHEM_DEBUG
             cout << Pe.ref->getIdentifier()<< "\t"<< Pe.comp->getIdentifier()<<endl;
#endif
           if (applyMOL2type)
           {
               const_cast<Atom*>(Pe.ref)->setMOL2Type(Pe.comp->getMOL2Type());

             for (vector<int>::const_iterator itMat2 =itMat+1; itMat2!= bestClic.end(); itMat2++)
             {
                 Pair &Pe2 = cliqueData.listPairs.at(*itMat2);
                 if (Pe.comp->hasBondWith(*Pe2.comp))
                 {
                    // if (verbose)  cout << PDB_name<<"|MATCHT|MODIFS\t"<<res->getIdentifier()<<"\tBOND\t"<< Pe.Ref->getAtom()->getName()<<"-"<<Pe2.Ref->getAtom()->getName()<<"\t"<<Pe.Comp->getAtom()->getBond(Pe2.Comp->getAtom())->getBondType()<<endl;
                     Bond *bde=const_cast<Bond*>(Pe.ref->getBondWith(*Pe2.ref));
                    const Bond *bdr=Pe.comp->getBondWith(*Pe2.comp);
                   // cout<< bdr<<endl;
                     bde->setBondType(bdr->getBondType());
                 }
             }
           }
             AtmDone.replace(residuOrder.at(Pe.ref),1,"1");
         }

         if (nAtChanged < residuOrder.size()) { cerr <<verbose<<"ERRORS\tNot the same atom counts "<< nAtChanged<<"\t"<<residuOrder.size()<<"\t"<<residu.getIdentifier()<<endl;}
         if ((residu.getResType() != ResType::STD_AA && templatedata.nHeavy > residuOrder.size())
            ||(templatedata.nHeavy > residuOrder.size()+1))
         {
             failedRes.push_back(&residu);
             cerr <<verbose<<"ERRORS\tINCOMPLETE RESIDU "<< templatedata.nHeavy<<"\t"<<residuOrder.size()<<"\t"<<residu.getIdentifier()<<endl;}

         if (AtmDone.find_first_of("0") != string::npos)
                 {
                     for (size_t n=0; n< AtmDone.length(); n++) if (AtmDone.substr(n,1)=="0") {

                         Atom* atmFailed=const_cast<Atom*>(atN.at(n));
                         atmFailed->setMOL2Type("Du");
                         failedAtom.push_back(atmFailed);
                         cerr  <<verbose<<"|ERRORS\t"<<atN.at(n)->getIdentifier()<<"\tunmatched"<<endl;
                     }
                 }

         for (ItCAtom itAA = residu.firstCAtom(); itAA != residu.lastCAtom();++itAA)
           {
             Atom &atm = **itAA;
            // cout << atm->getName()<< " " << atm->getName().substr(0,1)<<endl;
             if (atm.getName().substr(0,1)=="H")
               {

                 atm.setMOL2Type("H");
               }
           }

         if (residu.getResType() == ResType::STD_AA|| residu.getResType() == ResType::MOD_AA)
             for (ItCAtom itAA = residu.firstCAtom(); itAA != residu.lastCAtom();++itAA)
               {
                 Atom &atm = **itAA;
               // PATCH ERR::2
               const size_t NV = atm.getNumBond();
                 if (atm.getName() == "N")
                   {

                     if (NV == 2 || (NV==3 && atm.getResiduName()=="PRO") )                        atm.setMOL2Type("N.am");
                     else if (NV==3)
                     {
                        bool hasH=false;
                        for (size_t i=0;i<NV;++i) if (atm.getAtomLinked(i).isHydrogen())hasH=true;
                        if (hasH) atm.setMOL2Type("N.am");
                     }
                         else atm.setMOL2Type("N.4");
                   }
                 // END PATCH ERR::2
                 if (atm.getName() != "C" && atm.getName() != "SG") continue;

                 for (unsigned int nb=0;nb< NV;nb++)
                 {

                   Bond &bd= *atm.getBond(nb);
                   Atom &atmL=bd.getOtherAtom(atm);
                     if (atmL.getResidu() != &residu && atmL.getResidu()->getResType() == ResType::STD_AA &&atmL.getName() == "N")
                     {
                         bd.setBondType(BondType::AMIDE);
                         atmL.setMOL2Type("N.am");
                     }
                     else if  (atmL.getResidu() != &residu && atmL.getResidu()->getName() == "CYS" &&atmL.getName() == "SG")
                     {
                         bd.setBondType(BondType::SINGLE);
                     }
                     //PATCH ERR::2
                     else if (atmL.getName()=="O" && atm.getName()=="C")
                       {

                         bd.setBondType(BondType::DOUBLE);
                         atmL.setMOL2Type("O.2");

                       }
                     //END PATCH ERR::2
                 }
             }

        }catch (MoleExcept &e)
        {
            e.addTrace("Molecule::matchTemplate");
           e.addTrace("Involved residu: "+residu.getIdentifier());
           cerr << e.getCode()<<endl
                << e.getData()<<endl
                << e.getSource()<<endl
                << e.getTrace()<<endl;

        }
    }
}catch (MoleExcept &e)
{
        e.addTrace("Molecule::matchTemplate");
 cerr << e.getCode()<<endl
      << e.getData()<<endl
      << e.getSource()<<endl
      << e.getTrace()<<endl;
}
    delete[] templatesTable;

if (wVerbose)cout << verbose<<"END"<<endl;


}


void Molecule::cleanUnwanted(const std::string& verbose)
{

    const bool wVerbose= !verbose.empty();
    if (wVerbose)cout << verbose<<"START"<<endl;

  for (size_t i=0; i < Residues.size();++i)
  {
      Residu& residu= *Residues.at(i);

     // continue;
      if (((residu.getResType() == ResType::UNWANTED ||
           residu.getResType() == ResType::PROSTHETIC ||
           residu.getResType() == ResType::ION ||
           residu.getResType() == ResType::SUGAR) && residu.numAtom() != 0 && residu.getNInterResidu()==0) || residu.getResType() == ResType::ORGANOMET)
        {
          if (wVerbose)cout << verbose<<"MODIFS\t"<<residu.getIdentifier()<<"\t"<<residu.getLongName()<< "\t"<<residu.getNInterResidu()<<"\tDELETION\n";

          delResidu(residu);
     if (i>0)i--;
continue;
        }
      bool isFine=true;
      const size_t nAtm = residu.atoms.size();
      for (size_t j=0; j< nAtm;++j)
        {
          Atom& atom = *residu.atoms.at(j);
          try
          {
            atom.loadAtomicData();
          }
          catch (MoleExcept &e)
          {
              cerr << e.getCode()<<endl
                   << e.getData()<<endl
                   << e.getSource()<<endl
                   << e.getTrace()<<endl;

            isFine=false;
          }
        }
      if (!isFine)
        {
          if (wVerbose)cout << verbose<<"MODIFS\t"<<residu.getIdentifier()<<"\tUNWANTED ATOM\n";
          delResidu(residu);

        }
    }



if (wVerbose)cout << verbose<<"END\t"<<endl;

}







void Molecule::selChains(const std::string& ST,const std::string &verbose)
{
  vector<string> allowedchains;
  tokenStr(ST,allowedchains,"_");
  const bool wVerbose= !verbose.empty();
  if (wVerbose)cout << verbose<<"START"<<endl;

    ResiduList todel;
               double N_RES[NB_RESTYPE];
    const size_t nChain = Chains.size();
    for (size_t iCh=0; iCh < nChain;++iCh)
    {
        Chain& chain = *Chains.at(iCh);

        // Chain is temporary chain or allowed chain => skip it
        if (chain.getName() == "XX" ||
            find(allowedchains.begin(),allowedchains.end(),chain.getName()) != allowedchains.end())continue;
        // Check that chain is not a peptide
        for (size_t i=0; i < NB_RESTYPE;i++) N_RES[i]=0;
        const size_t nRes = chain.getCountRes();

        // Counting number of residue types within the chain :
        for(size_t iRes=0; iRes < nRes;++iRes)
        {
            N_RES[Residues.at(iRes)->getResType()]++;
        }
        // To be considered as peptide : less than 8 residues
        // with at least half of Standard or modified Amino Acid
        if (nRes <=8 && N_RES[ResType::STD_AA]+N_RES[ResType::MOD_AA] >= nRes/2)
        {
            continue;
        }
        if (wVerbose)cout << verbose<<"MODIFS\t"<<chain.getName()<<"\tDELETION\n";
        for(size_t iRes=0; iRes < nRes;++iRes)
        {
            todel.push_back(&chain.getResidu(iRes));
        }


    }// END SIZE_T iCh

    for (ItRes it = todel.begin(); it != todel.end(); it++)
        delResidu(**it);




if (wVerbose){ cout << verbose<<"END\t"<<endl;}
}


