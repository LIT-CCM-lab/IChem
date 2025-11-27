#include <iostream>
#include "headers/ICMole/atom.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/vertex.h"

using namespace ICMole;
using namespace std;


Residu::Residu(Molecule     &pmolecule,
               const  int         &Num,
               const  int         &FNum,
               const std::string  &Name,
               Chain              &chain,
               const unsigned int &ResType):
    identifier(""),
    name(Name),
    bits(""),
    num(Num),
    fileNum(FNum),
    resType(ResType),
    InterRes(0),
    rootAtom((Atom*)NULL),
    molecule(&pmolecule),
    chain(&chain),
    hetdata((HetData*)NULL),
    inUse(true)
{
    ostringstream oss;
    oss.width(3); oss << Name;
    oss.width(4);oss<<fileNum;
    oss<<"-";
    oss.width(2);oss<<chain.getName();
    chain.addResidu(this);
    identifier=oss.str();



}



Residu::Residu(Molecule *const    molecule,
               const  int         &Num,
               const  int         &FNum,
               const std::string  &Name,
               Chain              &chain,
               const unsigned int &ResType):
    identifier(Name),
    name(Name),
    bits(""),
    num(Num),
    fileNum(FNum),
    resType(ResType),
    InterRes(0),
    rootAtom((Atom*)NULL),
    molecule(molecule),
    chain(&chain),
    hetdata((HetData*)NULL),
    inUse(true)
{
    chain.addResidu(this);
    ostringstream oss;
    oss.width(3); oss << Name;
    oss.width(4);oss<<fileNum;
    oss<<"-";
    oss.width(1);oss<<chain.getName();
    identifier=oss.str();
}



/**
 * @brief Residu::addAtom
 * @param atom : atom to add
 *
 * Add an atom to this residu
 */
void Residu::addAtom(Atom& atom)
{
    atoms.push_back(&atom);
}

/**
 * @brief Residu::delAtom
 * @param atom : atom to remove from this residu
 *
 * Remove an atom to this residu
 */
void Residu::delAtom(const Atom& atom)
{
    const ItAtom it = std::find(atoms.begin(),atoms.end(), &atom);
    if (it != atoms.end()) atoms.erase(it);
}



/**
 * @brief Residu::getChainName
 * @return the name of the chain. Call chain->getName() function
 */
const std::string& Residu::getChainName() const
{
    return chain->getName();
}


/**
 * @brief Residu::toString
 * @return A string describing this residu, including atoms
 */
const std::string Residu::toString() const
{
    std::ostringstream oss;

    oss << "----------- RESIDU Num:"<< name <<"/"<<num<<"/"<<fileNum<<"/";
    if (chain != (Chain*)NULL) oss << chain->getName();
    oss << " ----------"<<endl;
    oss << " Molecule:"<< molecule->getName()<<" / "
        << " NINTERRES:"<< InterRes.size()<<std::endl
        << " ResType:"<< getLongName()<<" / "
        << " Bits : "<< bits<<endl;
    if (rootAtom != (Atom*)NULL) oss  << " ROOT :" <<rootAtom->getIdentifier()<<endl;


    if (atoms.size()){

        for (ItCAtom itA = atoms.begin(); itA != atoms.end(); itA++){

            oss << "    |-> " << (*itA)->getIdentifier()<<endl;

        }
    }
    return oss.str();

}

/**
 * @brief Residu::getLongName
 * @return A string describing the residu type (Amino Acid/Cofactor/Ligand ...)
 */
const std::string Residu::getLongName()const
{
    switch(resType)
    {
    case ResType::STD_AA:     return "Std Amino Acid";     break;
    case ResType::MOD_AA:     return "Modified Amino Acid";break;
    case ResType::NUCLEIC:    return "Nucleic  Acid";      break;
    case ResType::COFACTOR:   return "Cofactor";           break;
    case ResType::METAL:      return "Metal";              break;
    case ResType::ION:        return "Ion";                break;
    case ResType::WATER:      return "Water";              break;
    case ResType::LIGAND:     return "Ligand";             break;
    case ResType::PROSTHETIC: return "Prosthetic";         break;
    case ResType::ORGANOMET:  return "Organometallic";     break;
    case ResType::SUGAR:      return "Sugar";              break;
    case ResType::UNWANTED:   return "Unwanted";           break;

    }
    return "Unknown";
}




void Residu::checkUse()
{

    for (ItAtom itA= atoms.begin(); itA != atoms.end(); itA++)
    {
        if ((*itA)->isUsed()) {
            inUse=true;
            if (chain != (Chain*)NULL)chain->checkUse();
        }
    }
    inUse=false;
    if (chain != (Chain*)NULL)chain->checkUse();
}


/**
 * @brief Residu::setUse
 * @param newUse : True to keep this residu in the process. False otherwise
 * @param updateAtom : Update also its atom
 * @param toBonds : Update also atom bonds
 *
 */
void  Residu::setUse  (const bool& newUse,
                       const bool& updateAtom,
                       const bool& toBonds)
{
    inUse=newUse;
    if (!updateAtom)return;
    for (ItAtom itA= atoms.begin(); itA != atoms.end(); itA++)
    {
        (*itA)->setUse(newUse,false,toBonds);
    }
}






void Residu::checkResidu(bool wSetType) 
{








    if (!wSetType)return;
    for (unsigned int NAA = 0; NAA < NBAA; NAA++)
    {
        if (name == AAcid[NAA].name)
        {
            resType = AAcid[NAA].restype;

            break;
        }
    }



}




unsigned int Residu::Rules[NB_MOLETYPE][NB_RESTYPE];
bool Residu::rules_loaded=false;
void Residu::loadRules(const bool& force)
{
    if (rules_loaded && !force)return;
    rules_loaded=true;
    for (unsigned int I=0;I < NB_MOLETYPE;I++)
    {
        for (unsigned int J=0;J<NB_RESTYPE;J++)
        {
            Rules[I][J]=MoleType::UNDEFINED;
        }
    }
    Rules[MoleType::PROTEIN][ResType::STD_AA  ] = MoleType::PROTEIN;
    Rules[MoleType::PROTEIN][ResType::MOD_AA  ] = MoleType::PROTEIN;
    Rules[MoleType::PROTEIN][ResType::METAL   ]= MoleType::PROTEIN;
    
    Rules[MoleType::PROTEIN][ResType::WATER   ] = MoleType::UNDEFINED;
    Rules[MoleType::PROTEIN][ResType::COFACTOR] = MoleType::UNDEFINED;
    Rules[MoleType::PROTEIN][ResType::NUCLEIC ] = MoleType::PROTEIN;

    Rules[MoleType::LIGAND] [ResType::STD_AA  ] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::MOD_AA  ] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::SUGAR   ] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::ORGANOMET] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::UNWANTED] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::LIGAND  ] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::COFACTOR] = MoleType::LIGAND;
    Rules[MoleType::LIGAND] [ResType::UNKNOWN ] = MoleType::LIGAND;

    Rules[MoleType::COFACTOR] [ResType::COFACTOR ] = MoleType::LIGAND;
    Rules[MoleType::NUCLEIC][ResType::NUCLEIC ] = MoleType::UNDEFINED;

}


std::map<std::string,HetData> Residu::HETClass;
bool Residu::loadHETClass() 
{
    if (!HETClass.empty()) return true;
    try
    {
#ifdef ICHEM_DEBUG
        cout<< "ICHEM_DEBUG|"
            << "Loading HETClass file"
            <<endl;
#endif
        const std::string ICPath = get_IChem_LIB_Path();

        std::string ligne;
        std::vector<std::string> tokens;
        if (ICPath.length()==0)
            throw MoleExcept(1110101,
                             "GLOBAL::loadHETClass",
                             "No ICHEM_LIB DEFINED");


        std::ifstream ifs((ICPath+"/datas/HETLIST").c_str());
        if (!ifs.is_open())
            throw MoleExcept(1100102,
                             "GLOBAL::loadHETClass",
                             "Cannot open HETLIST file");

        while (!ifs.eof())
        {
            // GETTING DATA
            std::getline(ifs,ligne);
            tokens.clear();
            tokenStr(ligne,tokens,"\t");

            // ANALYSING DATA :
            if (tokens.size() != 5){continue;}
            HetData HD;
            HD.hasTemplate= (tokens.at(1) == "1") ? true:false;
            HD.HETreplace= (tokens.at(2) == "/") ? "": tokens.at(2);
            HD.MW = std::atof(tokens.at(3).c_str());
            const std::string& tmpType = tokens.at(4);


            if (tmpType == "STD_AA")         HD.ResT=ResType::STD_AA;
            else if (tmpType == "MOD_AA")         HD.ResT=ResType::MOD_AA;
            else if (tmpType == "COFACTOR")       HD.ResT=ResType::COFACTOR;
            else if (tmpType == "METAL")          HD.ResT=ResType::METAL;
            else if (tmpType == "WATER")          HD.ResT=ResType::WATER;
            else if (tmpType == "ION")            HD.ResT=ResType::ION;
            else if (tmpType == "LIGAND")         HD.ResT=ResType::LIGAND;
            else if (tmpType == "PROSTHETIC")     HD.ResT=ResType::PROSTHETIC;
            else if (tmpType == "ORGANOMET")      HD.ResT=ResType::ORGANOMET;
            else if (tmpType == "SUGAR")          HD.ResT=ResType::SUGAR;
            else if (tmpType == "UNWANTED")       HD.ResT=ResType::UNWANTED;
            else if (tmpType == "NUCLEIC")        HD.ResT=ResType::NUCLEIC;
            else                                  HD.ResT=ResType::UNKNOWN;
            HETClass.insert(std::pair<std::string,HetData>(tokens.at(0),HD));

        }
    }catch (MoleExcept &e)
    {
        if (e.getCode()==1110201) e.addTrace("GLOBAL::loadHETClass");
        throw;
    }
#ifdef ICHEM_DEBUG
    cout<< "ICHEM_DEBUG|End insertion"<<endl
        << "ICHEM_DEBUG|Number of HET Codes : "<< HETClass.size()
        <<endl;
#endif
    return true;
}



void Residu::applyResiduType(const std::string verbose) 
{
    const bool wVerbose= !verbose.empty();
    map<string,HetData>::iterator itHET, itHET2;

    try
    {
        if (Residu::HETClass.empty()) {
            Residu::loadHETClass();
        }
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Residu::applyResiduType");
        throw;
    };

    // Searching in the hetlist for the information
    itHET = HETClass.find(name);

    if (itHET == HETClass.end())
    {
        const string error = verbose+ "|ERROR\tNO DATA FOUND FOR : "+ identifier;
        switch (Moleaccess)
        {
        case Levels::NOTICE:cout<<error<<endl;return;break;
        case Levels::WARNING:cerr<<error<<endl;return;break;
        case Levels::STRICT:
        case Levels::FATAL:throw MoleExcept(3010101,"Residu::applyResiduType",error);
            break;
        }// END SWITCH
    }// END it END CLASS


    // Checking if name has been replace by a new one
    if ((*itHET).second.HETreplace != "")
    {
        if (wVerbose) cout << verbose <<"|REPLAC\t"  <<(*itHET).second.HETreplace<<"\n";
        // So we search for the new one
        itHET2 =HETClass.find((*itHET).second.HETreplace);

        // If the new one hasn't been found, then we keep the first one
        if (itHET2 == HETClass.end()){
            const string error =
                    verbose+
                    "|ERROR\tNO SUBSTITUTE "+
                    (*itHET).second.HETreplace+
                    " FOUND FOR : "+
                    identifier;
            switch (Moleaccess)
            {
            case Levels::NOTICE:cout<<error<<endl;break;
            case Levels::WARNING:cerr<<error<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:
                throw MoleExcept(3010102,
                                 "applyResiduType",
                                 error);
                break;
            }
            // Setting data for this residu
            hetdata=&(*itHET).second;
            // Setting the residu type
            resType=(*itHET).second.ResT;
        }

        // New HET found => apply this one
        else
        {
            hetdata=&(*itHET2).second;  //
            resType=(*itHET2).second.ResT;
        }


    }// END TEST HET REPLACE

    // No substitute name => apply the first one.
    else
    {
        hetdata=&(*itHET).second;
        resType=(*itHET).second.ResT;
        //  cout << "ICHEM_DEBUG: No replacement, ResT = " << resType << endl; // residue type = 5 here

        #ifdef ICHEM_DEBUG
            cout << "ICHEM_DEBUG|"<<verbose
            <<"|RESIDU\t"<<identifier
            << "\t"<<(*itHET).second.ResT
            << " " << getLongName()<<endl;
        #endif
    }// END ELSE


}



void Residu::addInterResidu( Residu* const residu)
{
    if (residu == (Residu*)NULL) return;
    InterRes.push_back(residu);
}


Atom& Residu::getAtom(const size_t&nAt) const
{
    return *atoms.at(nAt);
}



void Residu::getGraph(Graph &graph) const
{
    graph.clear();
    for (ItCAtom itA =atoms.begin();
         itA!=atoms.end();
         ++itA)

    {
        Atom& atom =**itA;
        graph.addVertex(atom.getNum(),atom.getName()).setAtom(&atom);

    }


    for (ItVert itV = graph.firstVertex();
         itV!= graph.lastVertex();
         ++itV)
    {
        Vertex &ve1 = **itV;
        for (ItVert itV2 = graph.firstVertex();
             itV2!= graph.lastVertex();
             ++itV2)
        {
            Vertex &ve2 = **itV2;
            if (ve1.getNum()>= ve2.getNum())continue;
            if (ve1.getAtom()->hasBondWith(*ve2.getAtom())) graph.addEdge(ve1,ve2);
        }
    }



}






void Residu::getDistMatrix(std::vector<unsigned short>& matrix,
                           std::map<const Atom*,unsigned short>& order,
                           const bool& wHydrogen) const
{
    try{
        std::map<const Atom*,unsigned short>::const_iterator itOrder,itOrderComp;

        matrix.clear();
        order.clear();
        unsigned short pos=0;
        for (ItCAtom itA =atoms.begin();
             itA!=atoms.end();
             ++itA)

        {
            const Atom& atom =**itA;

            if (!wHydrogen && atom.isHydrogen())continue;
            order.insert(pair<const Atom*,unsigned short>(&atom,pos));
#ifdef ICHEM_DEBUG
            cout << "TABLE ORDER : "<< atom.getName()<<"\t"<<pos<<endl;
#endif
            ++pos;

        }
        const size_t orderSize=order.size();
        matrix.reserve(orderSize*orderSize);
        // Filling it with dummy values (1000)
        matrix.assign(orderSize*orderSize,1000);

        // vertexdone is a boolean array used to avoid redundance during graph scan
        vector<bool>atomdone(orderSize,false);
        // vertextodo is the vector of vertex that we will currently look at
        // while vertextmp is the list of vertex that we will look at in the next loop
        vector<const Atom*> atomtodo, atomtmp;
        // Distance value from reference vertex to the other vertex
        unsigned short level=1;

        for (itOrder  = order.begin();
             itOrder != order.end();
             ++itOrder)
        {

            const Atom& atom = *(*itOrder).first;
#ifdef ICHEM_DEBUG
            cout << "REFERENCE ATOM : "<< atom.getName()<<endl;
#endif
            const unsigned short &posRef=(*itOrder).second;
            // Vertex against itself has a distance of 0
            matrix.at(posRef*orderSize+posRef)=0;

            // initialize test :
            atomtodo.clear();
            level=0;
            for (size_t pos=0; pos < orderSize;++pos) atomdone.at(pos)=false;
            atomdone.at(posRef)=true;

            // Scanning all linked vertex from reference vertex :
            for (size_t iRefVe=0;iRefVe < atom.getNumBond();++iRefVe)
            {
                const Atom &atomLink=atom.getAtomLinked(iRefVe);
#ifdef ICHEM_DEBUG
                cout << " |-> Linked atom : "<< atomLink.getName()<<endl;
#endif
                if (atomLink.getResidu() != this)continue;
                if (!wHydrogen && atomLink.isHydrogen())continue;

                itOrderComp=order.find(&atomLink);
                if (itOrderComp == order.end())
                    throw new MoleExcept(000000,
                                         "Residu::getDistMatrix",
                                         "Unexpected behavior");

                const unsigned short &posComp=(*itOrderComp).second;
#ifdef ICHEM_DEBUG
                cout << " |-> Position : "<< posComp<<endl;
#endif
                //if (order.at(std::distance(selectAtoms.begin(),)))
                // They have a distance of 1:
                matrix.at(posRef*orderSize+posComp)=1;
                matrix.at(posComp*orderSize+posRef)=1;
                // Pushing them into atomtodo so we can look at their edges
                atomtodo.push_back(&atomLink);
            }
#ifdef ICHEM_DEBUG
            cout << " |-> AtomToDo Size : "<< atomtodo.size()<<endl;
#endif

            while (!atomtodo.empty())
            {
                //  vertex to look in the next loop
                atomtmp.clear();

                level++;
#ifdef ICHEM_DEBUG
                cout << " |-> NEW LEVEL : "<< level<<endl;
#endif
                for (size_t iCompVe=0; iCompVe < atomtodo.size();++iCompVe)
                {
                    const Atom &atomComp = *atomtodo.at(iCompVe);
#ifdef ICHEM_DEBUG
                    cout << " |-> ATOM COMP : "<< atomComp.getName()<<endl;
#endif
                    if (!wHydrogen && atomComp.isHydrogen())continue;
                    itOrderComp=order.find(&atomComp);
                    if (itOrderComp == order.end())
                        throw new MoleExcept(000000,
                                             "Residu::getDistMatrix",
                                             "Unexpected behavior");
                    const unsigned short &posComp=(*itOrderComp).second;
#ifdef ICHEM_DEBUG
                    cout << " |-> Position COMP : "<< posComp<<endl;
#endif
                    // If we already did them (it's a closer path), we skip it
                    if (atomdone.at(posComp)) continue;
                    atomdone.at(posComp)=true;

                    // Setting up distance :
                    matrix.at(posRef*orderSize+posComp)=level;
                    matrix.at(posComp*orderSize+posRef)=level;
#ifdef ICHEM_DEBUG
                    cout << " |-> Distance : "<< atom.getName()<<"<->"<<atomComp.getName()<<"\t:"<<level<<endl;
#endif
                    // Fetching all edges of the comparison one:
                    for (size_t iCompVe=0;iCompVe < atomComp.getNumBond();++iCompVe)
                    {
                        const Atom &veLink=atomComp.getAtomLinked(iCompVe);

                        if (veLink.getResidu() != this)continue;
#ifdef ICHEM_DEBUG
                        cout << "     |-> Linked atom : "<< veLink.getName()<<" ";
#endif

                        itOrderComp=order.find(&veLink);

                        if (itOrderComp == order.end())
                        {
#ifdef ICHEM_DEBUG
                            cout<<endl;
#endif
                            continue;
                        }

#ifdef ICHEM_DEBUG
                        cout<<       atomdone.at((*itOrderComp).second)<<endl;
#endif

                        if (atomdone.at((*itOrderComp).second)) continue;

                        // And pushing them into atomtmp
                        atomtmp.push_back(&veLink);

                    }

                }
                // Switching atomtmp to atomtodo
                atomtodo.clear();
                if (atomtmp.empty())break;
                atomtodo=atomtmp;

            }// END WHILE


        }//END iRef

    }
    catch (MoleExcept &e)
    {
        cerr << e.getCode()<<endl
             << e.getData()<<endl
             << e.getSource()<<endl
             << e.getTrace()<<endl;
    }




}

const double& Residu::getWeight()
{
    const map<string,HetData>::const_iterator itHET = Residu::HETClass.find(name);
    if (itHET == Residu::HETClass.end())
    {
        throw MoleExcept(11111,"Complex::selectLigand","No HET Data found for this ligand");
    }
    return (*itHET).second.MW;
}
