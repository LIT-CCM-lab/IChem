#include <iostream>
#include <iomanip>
#include <ios>
#include "headers/ICPars/molereader.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/atom.h"
#include "headers/ICMole/bond.h"
//#define ICHEM_DEBUG // Debug pas refait, merci Jeremy
//# ICHEM_DEBUG2
using namespace std;
using namespace ICMole;


void MoleReader::loadNextMolecule(Molecule &molecule,
                                  const unsigned int& typeMole)
    try {

#ifdef ICHEM_DEBUG
//    if (debug)
        cout << string(80, '#') <<endl
             << setfill('#') << setw(80) << right << "NEW MOL2" << endl
             << string(80, '#') <<endl;
#endif

    molecule.setMoleType(typeMole);

    // CHECKS :
    if (!isOpen || !finput) throw MoleExcept(2010401,
                                             "MoleReader::loeadNextMolecule",
                                             "No file opened");
    if (finput.eof()) {EOFile=true;return;}

    int size_at   =0,   // Position in file of @<TRIPOS>ATOM
    size_bd   =0,   // Position in file of @<TRIPOS>BOND
    size_mo   =0,   // Position in file of the !!NEXT!! @<TRIPOS>MOLECULE
    size_set  =0,   // Position in file of @<TRIPOS>SET
    size_start=0;   // Beginning

    string ligne, name; // Each line of file goes here

    // Values set by @<TRIPOS>MOLECULE BLOCK
    unsigned int NbrAt   =100000,  // Number of atoms
            NbrBd   =0,       // Number of bonds
            NbrSubst=0,       // Number of substructure
            NbrFeat =0,       // Number of features
            NbrSets =0;       // Number of sets


    ////////////////////////////
    ////// READING HEADER //////
    ////////////////////////////
    /// Objective: knowing how many atom/bond/substructure we will have to
    /// handle
    size_start=finput.tellg();
    while ( !EOFile)
    {

        getLine(ligne);

        if (EOFile) break;
        if (ligne.find("@<TRIPOS>MOLECULE", 0)!=string::npos)
        {
            getLine(name);
            // GETTING THE NUMBER OF ATOMS/BONDS/SUBSTRUC :
            getLine(ligne);
            const int nscan = sscanf(ligne.c_str(),
                                     "%d %d %d %d %d",
                                     &NbrAt,
                                     &NbrBd,
                                     &NbrSubst,
                                     &NbrFeat,
                                     &NbrSets);
            if (nscan == EOF || nscan < 2)
                throw MoleExcept(2010402,
                                 "MoleReader::READ_MOL2",
                                 "Missing second data line");
            getLine(ligne);


            if (ligne.find("@<TRIPOS>") != string::npos)
                throw MoleExcept(2010403,
                                 "MoleReader::READ_MOL2",
                                 "Unexpected end of @<TRIPOS>MOLECULE block");

            // SETTING UP DATA :
            molecule.Atoms.reserve(NbrAt);
            molecule.Bonds.reserve(NbrBd);
            molecule.Residues.reserve(NbrSubst);
            molecule.name=name.substr(0,40);


#ifdef ICHEM_DEBUG
            cout << "MOLECULE" << endl;
            cout << setw(20) << left << "Molecule name"
                 << setw(10) << left << "N atoms"
                 << setw(11) << left << "N bonds"
                 << setw(11) << left << "N residues"
                 << setw(11) << left << "N features"
                 << setw(11) << left << "N sets"
                 << endl

                 << setw(20) << left << name
                 << setw(11) << left << NbrAt
                 << setw(11) << left << NbrBd
                 << setw(11) << left << NbrSubst
                 << setw(11) << left << NbrFeat
                 << setw(11) << left << NbrSets
                 << endl

                 << setw(15) << left << "Molecule type:" << typeMole << " (initial:"<<ligne<<")" << endl;
#endif
            // cout << "MOLECULE" << endl;
            // cout << setw(20) << left << "Molecule name"
            //     << setw(10) << left << "N atoms"
            //     << setw(11) << left << "N bonds"
            //     << setw(11) << left << "N residues"
            //     << setw(11) << left << "N features"
            //     << setw(11) << left << "N sets"
            //     << endl

            //     << setw(20) << left << name
            //     << setw(11) << left << NbrAt
            //     << setw(11) << left << NbrBd
            //     << setw(11) << left << NbrSubst
            //     << setw(11) << left << NbrFeat
            //     << setw(11) << left << NbrSets
            //     << endl

            //     << setw(15) << left << "Molecule type:" << typeMole << " (initial:"<<ligne<<")" << endl;
            break;
        } // END ligne - <TRIPOS>MOLECULE
    }// END WHILE !EOFILE


    ////////////////////////////
    ///////// CHECKING  ////////
    ////////////////////////////
    
    // Delete for the issue 
    // if (NbrAt==1000000) {
    //     EOFile=true;
    //     throw MoleExcept(2010404,
    //                      "MoleReader::READ_MOL2",
    //                      "NO MOLECULE FIND IN FILE :"+fileName);
    // }


    //////////////////////////////////////////////
    ///////// LOOKING FOR FILE POSITIONS  ////////
    //////////////////////////////////////////////
    /// To perform well, we need to first create substructures
    /// then atoms, bonds and finally sets
    /// So we need first to scan the whole file to fetch file positions
    /// and then analyse
    /// Since we need substructure first, we dont keep the position
    /// but directly read it
    while ( !EOFile )
    {
        getLine(ligne);
        if (finput.eof()) { EOFile=true;break; }
        if (ligne.find("@<TRIPOS>ATOM", 0)             !=string::npos)
        { size_at=finput.tellg();}
        if (ligne.find("@<TRIPOS>BOND", 0)             !=string::npos)
        { size_bd=finput.tellg();}
        if (ligne.find("@<TRIPOS>SET",  0)             !=string::npos)
        { size_set=finput.tellg();}
        if (ligne.find("@<TRIPOS>SUBSTRUCTURE", 0)     !=string::npos)
        { readMOL2Substruct(NbrSubst,molecule);}
        if (ligne.find("@<TRIPOS>MOLECULE")            != string::npos
                && size_at != 0 && size_bd != 0 )
        {finput.seekg(-ligne.length()-1,ios_base::cur);
            size_mo=finput.tellg();
            break;
        }
    }// END while EOFILE

    if (size_at==0)
        throw MoleExcept(2010405,
                         "MoleReader::READ_MOL2",
                         "No ATOM Block found");
    if (size_bd==0 && NbrBd!=0)
        throw MoleExcept(2010406,
                         "MoleReader::READ_MOL2",
                         "No BOND Block found");
    if (size_set==0 && NbrSets!=0)
        throw MoleExcept(2010407,
                         "MoleReader::READ_MOL2",
                         "No SET Block found");
    if (size_at == 0 && size_bd == 0)
    {
        EOFile=true;return;
    }
    if (finput.eof())  EOFile=true;


    ////////////////////////////
    ///////// ANALYSING  ///////
    ////////////////////////////
    /// Position file pointer to size_at to read atoms :
    finput.clear();
    finput.seekg(size_at,ios_base::beg);
    readMOL2Atom(NbrAt,molecule);
    getLine(ligne);
    if (ligne.find("@<TRIPOS>BOND") == string::npos)
        throw MoleExcept(2020408,
                         "MoleReader::READ_MOL2",
                         "Number of atoms count differs from actual atoms : "+ligne);

    /// Position file pointer to size_at to read bonds :
    finput.seekg(size_bd,ios_base::beg);
    if (NbrBd > 0)
        readMOL2Bond(NbrBd,molecule);
        
    getLine(ligne);
    if (ligne.find("@<TRIPOS>")==string::npos && ligne.find("#")==string::npos && ligne.length() > 0)
        throw MoleExcept(2020409,"MoleReader::READ_MOL2", "Premature end of bond block");

    finput.seekg(size_set,ios_base::beg);
    // if (NbrSubst > 0) readMOL2Set(NbrSets,molecule);
    //     if (NbrFeat > 0) readMOL2Feat(NbrFeat,size_start,size_mo,mole);
    if (size_mo != 0)
        finput.seekg(size_mo-1,ios_base::beg);

    ///////////////////////////////////
    ///////// POST-PROCESSING  ////////
    ///////////////////////////////////
    for (ItAtom itA = molecule.Atoms.begin(); itA != molecule.Atoms.end(); itA++)
    {
        try
        {

            (*itA)->checkMOL2type();
        }
        catch (MoleExcept &e)
        {
            switch (Moleaccess)
            {
            case Levels::NOTICE:cout<<e.getData()<<endl;break;
            case Levels::WARNING:cerr<<e.getData()<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:throw;break;
            }

        }
    }// END Scan Atoms

    if (typeMole == MoleType::LIGAND
            ||typeMole == MoleType::SITE
            ||typeMole == MoleType::COFACTOR
            ||typeMole == MoleType::WATER
            ||typeMole == MoleType::OTHER
            ||typeMole == MoleType::PROTEIN
            ||typeMole == MoleType::UNDEFINED
            ||typeMole == MoleType::NUCLEIC)
    {
        bool HETLoaded=false;
        try
        {
            HETLoaded=Residu::loadHETClass();
        }catch (...){};
        try{


            for (ItRes res =molecule.Residues.begin();
                 res!=molecule.Residues.end();
                 ++res)
            {
                if (HETLoaded) (*res)->applyResiduType();
                (*res)->checkResidu(!HETLoaded);
                //            cout << (*res)->getName() << endl;
            }

        }
        catch (MoleExcept &e)
        {
            e.addTrace("MoleReader::READ_MOL2");
            throw;
        }
    }



}
catch (MoleExcept &e)
{
    e.addTrace("MoleReader::loadNextMolecule");
    throw;
}









/*! \fn void MoleReader::Read_MOL2_Substruct(const unsigned int NbrSubst, Molecule* mole)
* \param NbrSubst : Number of substructure to search
* \param mole : Molecule to put atom in
*
* \todo Handling when no  subst_subType is empty
*/
void MoleReader::readMOL2Substruct(const unsigned int NbrSubst,
                                   Molecule& mole)
throw(MoleExcept)
{
#ifdef ICHEM_DEBUG
    cout << "__________________________________________________"<<endl
         << "------------------ SUBSTRUCTURES ------------------"<<endl;
#endif

    string ligne;


    string subst_Name   ="", // Name of the residue i.e. ASP168
            subst_Type   ="", // Type of the residue TEMP, PERM, RESIDUE, GROUP or DOMAIN
            subst_Chain  ="", // Chain of the residue i.e. A
            subst_subType="", // SubType of the residu i.e. ASP
            subst_Bits   =""; // Sybyl BITS
    unsigned int subst_AtmR   =0,  // Id of the root atom i.e. 1042
            subst_IntBd  =0,  // Number of inter substructure bonds
            subst_Pos    =0,  // If of the residue, given by file
            selId        =0,  // Id of the residu, took by deleting subType to Name i.e. 168
            np           =0;  // Used for splitting line
    ostringstream oss;              // Used for mapping
    mapRes.clear();
    mapRes.insert(pair<string,Residu*>("100000TMP",&mole.getUnknownRes()));

    for (unsigned int NSub=1; NSub <= NbrSubst; NSub++)
    {
        getLine(ligne);
        // Checking if it's a comment => passing
        if (ligne.substr(0,1) == "#") {NSub--;continue;}
#ifdef ICHEM_DEBUG
        cout << "ICHEM_DEBUG|SUBSTRUCT --"<<ligne<<endl;
#endif
        if (ligne == "END" || ligne.find("@<TRIPOS>", 0)!=string::npos)
            throw MoleExcept(2020501,
                             "MoleReader::readMOL2Substruct",
                             "Unexpected end of substucture block");

        // Cleaning data :
        subst_Name="", subst_Type="",subst_Chain="";subst_Bits="";
        subst_AtmR=0; subst_IntBd=0;subst_subType="";selId=0;

        // Skip delimiters at beginning.
        string::size_type lastPos = ligne.find_first_not_of(" \t", 0);
        // Find first "non-delimiter".
        string::size_type pos     = ligne.find_first_of(" \t", lastPos);
        np=0;
        while (string::npos != pos || string::npos != lastPos)
        {
            np++;
            // Found a token, add it to the corresponding position
            switch (np)
            {
            case 1: subst_Pos     = atoi(ligne.substr(lastPos,pos-lastPos).c_str());break;
            case 2: subst_Name    =      ligne.substr(lastPos,pos-lastPos);      break;
            case 3: subst_AtmR    = atoi(ligne.substr(lastPos,pos-lastPos).c_str());break;
            case 4: subst_Type    =      ligne.substr(lastPos,pos-lastPos);         break;
                // case 5 => Type of directory - unnecessary
            case 6: subst_Chain   =      ligne.substr(lastPos,pos-lastPos);         break;
            case 7: subst_subType =      ligne.substr(lastPos,pos-lastPos);         break;
            case 8: subst_IntBd   = atoi(ligne.substr(lastPos,pos-lastPos).c_str());break;
            case 9: subst_Bits    =      ligne.substr(lastPos,pos-lastPos);         break;
            }

            // Skip delimiters.  Note the "not_of"
            lastPos = ligne.find_first_not_of(" \t", pos);
            // Find next "non-delimiter"
            pos = ligne.find_first_of(" \t", lastPos);
        }
        if (np < 3) throw MoleExcept(2020502,
                                     "MoleReader::readMOL2Substruct",
                                     "Not enough descriptors for substructure : "
                                     +ligne);



#ifdef ICHEM_DEBUG
        cout << "ICHEM_DEBUG|" << "Name:"<<subst_Name<<"  \tRootAtm:"<<subst_AtmR<<"  \tType:"<<subst_Type<<endl
             << "Chain:"<<subst_Chain<<"  \tSubType:"<<subst_subType<<"  \tIntDb:"<<subst_IntBd<<endl;
#endif

        try
        {
            // No Given name - Atoms goes to unknownResidu:
            if (subst_Name == "****" ||subst_Name=="") subst_Name="TMP";
            // No given chain - Residues goes to unknownChain
            if (subst_Chain == "****"||subst_Chain=="") subst_Chain="XX";
            // Otherwise :
            if (subst_subType != "" && subst_subType != "****")
            {

                if (subst_Name.length() > subst_subType.length())
                {
                    selId = atoi(subst_Name.substr(subst_subType.length()).c_str());
                    if (subst_Name.find(subst_subType)==string::npos)
                    {
                        ostringstream error;
                        error << "subst_name differ from sub_type for"<<endl
                              << ligne<<endl;
                        switch (Moleaccess)
                        {
                        case Levels::NOTICE:cout <<"WARNING|"<<error.str();break;
                        case Levels::WARNING:cerr <<"WARNING|"<<error.str();break;
                        case Levels::STRICT:
                        case Levels::FATAL:
                            throw MoleExcept(2020503,
                                             "MoleReader::readMOL2Substruct",
                                             error.str());
                        }
                    }
                }
                else
                {

                    // Handle nucleic acid :
                    if ((subst_Name.substr(0,1) == "G" ||
                         subst_Name.substr(0,1) == "C" ||
                         subst_Name.substr(0,1) == "A" ||
                         subst_Name.substr(0,1) == "T")
                            && isdigit(subst_Name[1])
                            )
                        subst_subType=subst_Name[0];
                    selId = atoi(subst_Name.substr(1).c_str());
                }
            }
            else { // When issues with name
                // We scan AACid list to fetch for any similarity
                bool find_RT=false;
                for (unsigned int NAA = 0; NAA < NBAA; NAA++)
                {

                    if (subst_Name.find(AAcid[NAA].name) != string::npos)
                    {
                        find_RT=true;
                        if (subst_Name.length() > AAcid[NAA].name.length())
                            selId = atoi(subst_Name.substr(AAcid[NAA].name.length()).c_str());
                        subst_subType=AAcid[NAA].name;
                    }
                }

                if (!find_RT)
                {
                    if (subst_Name.length() > 3) subst_subType=subst_Name.substr(0,3);
                    else subst_subType=subst_Name;
                }
            }

            Residu &res = mole.addResidu(subst_Chain,selId,subst_subType);

            if (subst_AtmR != 0) mapResRoot.insert(pair<unsigned int,Residu*>(subst_AtmR,&res));
            res.InterRes.reserve(subst_IntBd);
            res.bits= subst_Bits;
            oss.str(""); oss << subst_Pos<<subst_Name;
            mapRes.insert(pair<std::string,Residu*>(oss.str(),&res));
            //         cout << subst_Pos<<subst_Name<< endl;
        }
        catch (MoleExcept &e)
        {
            if (e.getCode() != 2020503
                    &&e.getCode() != 2020502
                    &&e.getCode() != 2020501)
                e.addTrace("MoleReader::readMOL2Substruct");
            throw;
        }

    }// END FOR NSub


}


void  MoleReader::readMOL2Atom(const unsigned int NbrAt,
                               Molecule &mole) throw(MoleExcept)
{
#ifdef ICHEM_DEBUG
    cout << "ICHEM_DEBUG|____________________________________________"<<endl
         << "ICHEM_DEBUG|------------------- ATOMS ------------------"<<endl;
#endif

    string ligne;

    unsigned int atm_resId   = 0;
    unsigned int atm_Id      = 0;
    string       atm_resName = "";
    string       atm_Name    = "";

    std::map<std::string,Residu*>::iterator ITTI;
    ostringstream oss;
    Residu *res;
    std::map<unsigned int,Residu*>::iterator ITTR;
    unsigned int np=0;
    string::size_type lastPos;
    string::size_type pos;
    Atom* atom=(Atom*)NULL;
    mole.bary_check=false;
    for (unsigned int NAt = 1; NAt <= NbrAt; NAt++)
    {
        getLine(ligne);
        if (ligne.substr(0,1) == "#") {NAt--;continue;}
#ifdef ICHEM_DEBUG
        cout << "ICHEM_DEBUG| -- ATOM --"<<ligne<<endl;
#endif
        if (ligne == "END" || ligne.find("@<TRIPOS>", 0) != string::npos)
        {
            oss.str("");
            oss<< "Unexpected end of atom block. Expecting "
               << NbrAt
               << " atoms. Found "
               << NAt-1<<" atoms"<<endl;
            throw MoleExcept(2020601,"MoleReader::readMOL2Atom",oss.str());
        }

        atom=(Atom*)NULL;

        atm_resName="";
        atm_resId=0;
        try{

            atom = new Atom();
            atom->molecule=&mole;




            // Skip delimiters at beginning.
            lastPos = ligne.find_first_not_of(" \t", 0);
            // Find first "non-delimiter".
            pos     = ligne.find_first_of(" \t", lastPos);
            np=0;
            while (string::npos != pos || string::npos != lastPos)
            {
                np++;
                // Found a token, add it to the vector.
                switch (np)
                {
                case 1: atom->fNum  = atoi(ligne.substr(lastPos,pos-lastPos).c_str());break;
                case 2: atom->name    = ligne.substr(lastPos,pos-lastPos);break;
                case 3: atom->fixpos.x       = round(1000*atof(ligne.substr(lastPos,pos-lastPos).c_str()))/1000;break;
                case 4: atom->fixpos.y       = round(1000*atof(ligne.substr(lastPos,pos-lastPos).c_str()))/1000;break;
                case 5: atom->fixpos.z       = round(1000*atof(ligne.substr(lastPos,pos-lastPos).c_str()))/1000;break;
                case 6: atom->mol2type    = ligne.substr(lastPos,pos-lastPos);break;
                case 7: atm_resId   = atoi(ligne.substr(lastPos,pos-lastPos).c_str());break;
                case 8: atm_resName = ligne.substr(lastPos,pos-lastPos);break;
                case 9: atom->partial_charge  = atof(ligne.substr(lastPos,pos-lastPos).c_str());break;
                case 10:atom->bits        = ligne.substr(lastPos,pos-lastPos);break;
                }

                


                
                
                 
                
                


                // Skip delimiters.  Note the "not_of"
                lastPos = ligne.find_first_not_of(" \t", pos);
                // Find next "non-delimiter"
                pos = ligne.find_first_of(" \t", lastPos);
            }

            // cout << "id_atom:  " << atom->fNum << " | " << "atom_name: " << atom->name  << " | " << "x: " << atom->fixpos.x  << " | " << "y: " << atom->fixpos.y << " | " << 
            //     "z: " << atom->fixpos.z << " | " << "mol2type:  " << atom->mol2type << " | " << "atom_res_id: " << atm_resId << " | " << "atom_res_name: " << atm_resName << " | " <<
            //     "partial_charge: " << atom->partial_charge << endl;

            if (np < 6)
                throw MoleExcept(2020602,
                                 "MoleReader::readMOL2Atom",
                                 "Not enough descriptors for atom : "+ligne);





        }
        catch (std::bad_alloc& ba)
        {

            Moleaccess=Levels::FATAL;
            cerr << ba.what()<<endl;
            throw MoleExcept(1070102,
                             "Molecule::addAtom",
                             "Bad allocation exception - exiting");
        }
        catch (MoleExcept &e)
        {
            delete atom;
            throw;
        }
        mole.Atoms.push_back(atom);
        atom->num=mole.maxNumAtom;
        mole.maxNumAtom++;




#ifdef ICHEM_DEBUG
        cout << "ICHEM_DEBUG|ID:"       <<atom->fNum
             <<"  \tName:"  <<atom->name
            <<"  \tMOL2:"  <<atom->mol2type
           <<"  \tCharge:"<<atom->partial_charge<<endl
          << "ICHEM_DEBUG|X:"<<atom->fixpos.x
          <<"  \tY:"<<atom->fixpos.y
         <<"  \tZ:"<<atom->fixpos.z<<endl
        << "ICHEM_DEBUG|ResID:"<<atm_resId
        <<"  \tResName:"<<atm_resName<<endl
        <<"ICHEM_DEBUG|#########"<<endl;
#endif

        // Analysing value:
        const std::string former_name=atm_resName;
        if (atm_resName.substr(0,1) == "<" && atm_resName.substr(atm_resName.length()-1) == ">") {
            atm_resName="TMP"; 
            atm_resId=1;
        }
        else if (atm_resId == 0 && atm_resName=="") {
            atm_resName = "TMP";
            atm_resId = 100000;
        }

        //Searching corresponding residu :
        oss.str("");
        oss << atm_resId << atm_resName;
        ITTI = mapRes.find(oss.str());

        if (ITTI == mapRes.end()) {
            ostringstream error;
            cerr << "Residue " << former_name << "/" << atm_resId << " was not found for the atom " << atm_Name << "/" << atm_Id;
            
            switch (Moleaccess) {
            case Levels::NONE:
                break;
            case Levels::NOTICE:
                cerr << "WARNING|" << error.str() << "\n"; 
                break;
            case Levels::WARNING:
                cerr << "ERROR|" << error.str() << "\n";
                break;
            case Levels::STRICT:
            case Levels::FATAL:
                throw MoleExcept(2020703,"MoleReader::readMOL2Atom", error.str());
                break;
            }// END SWITCH

            res = &mole.unknownRes;
        }
        else res=(*ITTI).second;

        atom->residu=res;
        res->atoms.push_back(atom);


        ITTR = mapResRoot.find(atom->fNum);
        if (ITTR != mapResRoot.end()
                && res != &mole.unknownRes
                && atom->getName()=="CA")
        {

            (*ITTI).second->rootAtom=atom;
        }

    }


}



/*! \fn void MoleReader::Read_MOL2_Bond(const unsigned int NbrBd, Molecule* mole) throw(MoleExcept)
* \param NbrBd : Number of bonds to search
* \param mole : Molecule to put atom in
*
*/
void MoleReader::readMOL2Bond(const unsigned int NbrBd, Molecule& mole) throw(MoleExcept)
{
#ifdef ICHEM_DEBUG
    cout << "ICHEM_DEBUG|____________________________________________"<<endl
         << "ICHEM_DEBUG|------------------- BONDS ------------------"<<endl;
#endif


    const unsigned int sizeAtm = (unsigned int)mole.Atoms.size();
    unsigned int idAtm1 =0, // First atom involved in the bond
            idAtm2 =0, // Second atom involved in the bond
            bdId   =0, // Bond ID
            np     =0;
    string       idBType="",// Bond Types
            bits   ="",// Sybyl Bits
            ligne  ="";
    short        BondTypes;

    for (unsigned int NBd=1; NBd <= NbrBd; NBd++)
    {
        getLine(ligne);
        if (ligne.substr(0,1) == "#") {NBd--;continue;}
#ifdef ICHEM_DEBUG
        cout << "ICHEM_DEBUG| -- BOND --"<<ligne<<endl;
#endif
        if (ligne == "END" || ligne.find("@<TRIPOS>", 0)!=string::npos)
            throw MoleExcept(2020701,
                             "MoleReader::readMOL2Bond",
                             "Unexpected end of bond block");

        // Reinitialize values :
        bdId    =0;
        idAtm1  =0;
        idAtm2  =0;
        idBType="";
        bits="";
        // Skip delimiters at beginning.
        string::size_type lastPos = ligne.find_first_not_of(" \t", 0);
        // Find first "non-delimiter".
        string::size_type pos     = ligne.find_first_of(" \t", lastPos);
        np=0;
        while (string::npos != pos || string::npos != lastPos)
        {
            np++;
            // Found a token, add it to the vector.
            switch (np)
            {
            case 1 : bdId = atoi(ligne.substr(lastPos,pos-lastPos).c_str());break;
            case 2: idAtm1 = atoi(ligne.substr(lastPos,pos-lastPos).c_str())-1;break;
            case 3: idAtm2 = atoi(ligne.substr(lastPos,pos-lastPos).c_str())-1;break;
            case 4: idBType = ligne.substr(lastPos,pos-lastPos);break;
            case 5: bits=ligne.substr(lastPos,pos-lastPos);break;
            }

            // Skip delimiters.  Note the "not_of"
            lastPos = ligne.find_first_not_of(" \t", pos);
            // Find next "non-delimiter"
            pos = ligne.find_first_of(" \t", lastPos);
        }

        // if (np < 4)
        //     throw MoleExcept(2020702,
        //                      "MoleReader::readMOL2Bond",
        //                      "Not enough parameter in bond block : "+ligne);
        // if (idAtm1 >= sizeAtm )
        //     throw MoleExcept(2020703,
        //                      "MoleReader::readMOL2Bond",
        //                      "Atom 1 is above atom number on line :"+ligne);
        // if (idAtm2 >= sizeAtm )
        //     throw MoleExcept(2020704,
        //                      "MoleReader::readMOL2Bond",
        //                      "Atom 2 is above atom number on line :"+ligne);

        Atom *const atm1 = mole.Atoms.at(idAtm1);
        Atom *const atm2 = mole.Atoms.at(idAtm2);


        if (idBType.find("1")!=string::npos)       BondTypes = BondType::SINGLE;
        else if (idBType.find("2")!=string::npos)  BondTypes=  BondType::DOUBLE;
        else if (idBType.find("3")!=string::npos)  BondTypes=  BondType::TRIPLE;
        else if (idBType.find("ar")!=string::npos) BondTypes=  BondType::AROMATIC;
        else if (idBType.find("am")!=string::npos) BondTypes = BondType::AMIDE;
        else
        {
            switch (Moleaccess)
            {
            case Levels::NOTICE:
                cout << atm1->getIdentifier()<<"\t"<<
                        atm2->getIdentifier()<<"\tDUMMY BOND"<<endl;break;
                break;
            case Levels::WARNING:
                cerr << atm1->getIdentifier()<<"\t"<<
                        atm2->getIdentifier()<<"\tDUMMY BOND"<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:
                throw MoleExcept(2020705,
                                 "MoleReader::readMOL2Bond",
                                 atm1->getIdentifier()+"\t"+
                                 atm2->getIdentifier()+"\tDUMMY BOND");
                break;
            }
            BondTypes = BondType::DUMMY;
        }
        // cout << "show me the bonds " << BondTypes << endl;
        Bond &bd=mole.addBond(atm1,atm2,BondTypes);
        bd.setBits(bits);

    }



}
