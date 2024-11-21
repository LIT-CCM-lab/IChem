#include <iostream>
#include "headers/ICPars/molewriter.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICMole/bond.h"
#include "headers/ICMole/complex.h"
using namespace std;
using namespace ICMole;

MoleWriter::MoleWriter():
    fname(""),
    FFormat(FileFormat::MOL2),
    withH(true),
    withRotPos(false),
    withCycles(false),
    withSets(true),
    onlyUsed(false)

{}


MoleWriter::MoleWriter(const string &name,
                       const short &FFormat):
    fname(name),
    FFormat(FFormat),
    withH(true),
    withRotPos(false),
    withCycles(false),
    withSets(true),
    onlyUsed(false)
{

    if (fname != "")  openFile();
}



/**
 * @brief MoleWriter::newFile
 * @param newfname: path of the file
 * @param FForma : File format
 * @return True if the file has been open. False otherwise
 *
 * Close a potentially already opened file, and open the new one defined by fname
 */
bool MoleWriter::newFile(const string& newfname, const unsigned int& FForma)
{

    if (ofstr.is_open())ofstr.close();
    fname=newfname;FFormat=FForma;
    return openFile();
}


/**
 * @brief MoleWriter::openFile
 * @return  True if the file has been open. False otherwise
 *
 * Open the file
 */
bool MoleWriter::openFile()
{

    if (fname.length() ==0)return false;
    ofstr.open(fname.c_str(),ios::out|ios::binary);
    if (!ofstr.is_open())
    {    cout <<"Error Writing to "<< fname <<endl;return false;}
    return true;
}

MoleWriter::~MoleWriter() { if (ofstr.is_open()) ofstr.close();}


void MoleWriter::printMOL2(Molecule *const mole) {
    if (mole== (Molecule*)NULL)return;

    time_t timestamp = time(NULL);
    char buffer[256];
    unsigned int Natm=0;
    unsigned int  Nbd=0;
    unsigned int NRes=0;
    unsigned int NSet=0;
    unsigned int  RId=0;
    int RNum=0;
    unsigned int AccVal=0;
    string RName="";

    const size_t maxAtm= mole->getMaxNumAtom();
    const size_t maxCyclAtom=mole->numCyclAtom();
    const size_t maxAromAtom=mole->numAromAtom();
    const size_t maxUnkAtom =mole->numUnkAtom();
    const size_t nAllAtom = maxAtm+maxCyclAtom+maxAromAtom+maxUnkAtom;

    size_t AtmAccess[nAllAtom];
    size_t BdAccess[mole->getMaxNumBond()];
    size_t ResAccess[mole->getMaxNumResidu()+3];
    Coords coo;


    mole->getCycleRes().setNum(mole->getMaxNumResidu()+1);
    mole->getAromaticRes().setNum(mole->getMaxNumResidu()+2);
    mole->getUnknownRes().setNum(mole->getMaxNumResidu()+3);


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// INITIALIZATION ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    for (size_t i=0;i<nAllAtom;i++){AtmAccess[i]=0;   }
    for (size_t i=0;i<mole->getMaxNumBond();i++){BdAccess[i]=0;}
    for (size_t i=0;i<mole->numResidus()+3;i++){ResAccess[i]=0;}


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////// PRE-CALCULATIONS /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;

        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {

            Atom *atmNP = *itA;

            if (!withH && atmNP->isHydrogen()){AtmAccess[atmNP->getNum()]=0;continue; }
            if (onlyUsed && !atmNP->isUsed()){AtmAccess[atmNP->getNum()]=0;continue; }
            Natm++;
            //cout << atmNP->getIdentifier()<<" " <<atmNP->getNum()<<" " << Natm<<endl;
            AtmAccess[atmNP->getNum()]=Natm;
        }
    }

    for (ItCBond itB = mole->firstBond();
         itB != mole->lastBond();
         itB++)
    {
        try{
            if ((AtmAccess[(*itB)->getAtom1().getNum()] == 0
                 ||AtmAccess[(*itB)->getAtom2().getNum()] == 0)
                || (onlyUsed && !(*itB)->isUsed()))
            {
                BdAccess[(*itB)->getNum()]=0;
                continue;
            }
        }
        catch (...)
        {
            cerr << "ERRS : "<< (*itB)->getAtom1().getIdentifier()<<"\t"<<(*itB)->getAtom2().getIdentifier()<<endl;
        }
        Nbd++;
        BdAccess[(*itB)->getNum()]=Nbd;
    }
    bool auth=false;
    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {
        auth=false;
        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {

            if (AtmAccess[(*itA)->getNum()] != 0){auth=true;break;}
        }
        if (!auth) {ResAccess[res->getNum()]=0;continue;}
        NRes++;
        ResAccess[res->getNum()]=NRes;


    }
    if (maxUnkAtom){NRes++;ResAccess[mole->getCycleRes().getNum()]=NRes;}
    if (withCycles)
    {
        if (maxCyclAtom >0){NRes++;ResAccess[mole->getCycleRes().getNum()]=NRes;}
        if (maxAromAtom >0){NRes++;ResAccess[mole->getAromaticRes().getNum()]=NRes;}
        NSet+=maxCyclAtom+maxAromAtom;
    }

    ostringstream ossCol;
    if (withSets)  { NSet+=mole->getNumSet();

        vector<unsigned int> TablCol[NB_COLOR];
        for(ItCAtom itAS = mole->firstColAtom(); itAS != mole->lastColAtom(); itAS++)
        {
            const Atom& atm = **itAS;
            if (AtmAccess[atm.getNum()]==0)continue;


            for (unsigned int ICol = 0; ICol < NB_COLOR;ICol++)
            {

                if (atm.getRcolor() == ColorData[ICol].RED
                    && atm.getGcolor() == ColorData[ICol].GREEN
                    && atm.getBcolor() == ColorData[ICol].BLUE)
                {

                    TablCol[ICol].push_back(AtmAccess[atm.getNum()]);break;}
            }
        }

        for (unsigned int ICol=0; ICol < NB_COLOR; ICol++)
        {

            if (TablCol[ICol].size()== 0) continue;
            NSet++;
            ossCol<< "ATOM$"<<ColorData[ICol].name<<" STATIC     ATOMS   COLORGROUP SYSTEM "<<endl;
            ossCol<< TablCol[ICol].size()<<" ";
            unsigned int count=1,countTOT=0;
            for (vector<unsigned int>::iterator it = TablCol[ICol].begin(); it != TablCol[ICol].end(); it++)
            {
                ossCol<< *it<<" ";count++;countTOT++;
                if (count == 18&& countTOT != TablCol[ICol].size()) {count=0; ossCol <<"\\"<<endl;}
            }
            ossCol<<endl;
        }


    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// HEADER BLOCK ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////


    cout<< "#    Name: "<< mole->getName()<<"\n";
    strftime(buffer, sizeof(buffer), "%c", localtime(&timestamp));
    cout << "#    Modification time : " << buffer<<"\n"<<"\n";
    cout <<"@<TRIPOS>MOLECULE"<<"\n";
    cout<< mole->getName()<<"\n";
    cout.setf(ios::right);cout.width(5); cout<<Natm<<" ";
    cout.setf(ios::right);cout.width(5); cout<<Nbd<<" ";
    cout.setf(ios::right);cout.width(5); cout<<NRes<<" ";
    cout.setf(ios::right);cout.width(5); cout<<"0"<<" ";
    cout.setf(ios::right);cout.width(5); cout<<NSet<<"\n";

    // Molecular name :
    switch (mole->getMoleType())
    {
        case MoleType::PROTEIN: cout<< "PROTEIN"<<"\n";break;
        case MoleType::COFACTOR:
        case MoleType::LIGAND:
        case MoleType::WATER:
        default: cout<<"SMALL"<<"\n";break;
    }

    cout << "NO_CHARGES"<<"\n"<<"\n"<<"\n";
    if (mole->getMoleType() == MoleType::PROTEIN)
        cout<< "@<TRIPOS>DICT"<<"\n"<<"PROTEIN macromol"<<"\n";


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// ATOM BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    cout<< "@<TRIPOS>ATOM"<<"\n";
    cout.setf(ios_base::fixed,ios_base::floatfield);cout.precision(4);
    unsigned int IdAtom,IdBond;
    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {
            Atom *atm=*itA;
            IdAtom=AtmAccess[atm->getNum()];
            if (IdAtom==0 )continue;

            RName=atm->getResidu()->getName();
            RId=ResAccess[atm->getResidu()->getNum()];
            RNum=atm->getResidu()->getFNum();


            if (withRotPos) coo = atm->rotpos; else coo = atm->fixpos;

            cout.setf(ios::right);cout.width(7);cout<<IdAtom<<" ";
            cout.unsetf(ios::right);cout.setf(ios::left); cout.width(7);cout<<atm->getName()<<" ";
            cout.unsetf(ios::left); cout.setf(ios::right);
            cout.width(10);cout<<coo.x<<" ";
            cout.width(9); cout<<coo.y<<" ";
            cout.width(9); cout<<coo.z<<" ";
            cout.width(6); cout.unsetf(ios::right);cout.setf(ios::left);  cout<<atm->getMOL2Type()<<" ";
            cout.width(4); cout.unsetf(ios::left); cout.setf(ios::right); cout<<RId<<" ";
            if (RName.length() >=3)
            {
                cout.width(3); cout.unsetf(ios::right);cout.setf(ios::left);   cout<<RName;
                cout.width(5);                                                  cout<<RNum<<" ";
            }
            else if (RName.length()==2)
            {
                cout.width(2); cout.unsetf(ios::right);cout.setf(ios::left);   cout<<RName;
                cout.width(6);                                                  cout<<RNum<<" ";

            }
            cout.width(9); cout.unsetf(ios::left); cout.setf(ios::right); cout<<atm->getPartialCharge()<<" ";
            if (atm->getBits().length()) cout << atm->getBits();
            cout<<"\n";
        }}

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    cout<< "@<TRIPOS>BOND"<<"\n";
    for (ItCBond it = mole->firstBond(); it != mole->lastBond();it++)
    {
        Bond *bd=*it;
        IdBond=BdAccess[bd->getNum()];
        if (IdBond == 0) continue;
        cout.setf(ios::right);cout.width(6);cout<<IdBond<<" ";
        cout.width(4); cout<<AtmAccess[bd->getAtom1().getNum()]<<" ";
        cout.width(4); cout<<AtmAccess[bd->getAtom2().getNum()]<<" ";
        cout.unsetf(ios::right);cout.setf(ios::left);cout.width(4);
        switch (bd->getBondType())
        {
            case BondType::SINGLE: cout<<"1";break;
            case BondType::DOUBLE: cout<<"2";break;
            case BondType::TRIPLE: cout<<"3";break;
            case BondType::AMIDE:  cout<<"am";break;
            case BondType::AROMATIC:  cout<<"ar";break;
            case BondType::ANY:  cout<<"an";break;
            case BondType::DUMMY:  cout<<"du";break;
            default:  cout<<"un";break;
        }
        cout<<" ";
        if (bd->getBits().length()) cout << bd->getBits();
        cout<<"\n";
    }


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    string currChain="";

    cout<< "@<TRIPOS>SUBSTRUCTURE"<<"\n";
    ostringstream oss;
    for (ItCRes it = mole->firstResidu(); it != mole->lastResidu();it++)
    {
        Residu *res = *it;
        AccVal=ResAccess[res->getNum()];
        if (AccVal == 0) continue;
        cout.setf(ios::right);cout.width(6);cout<<AccVal<<" ";
        cout.unsetf(ios::right);cout.setf(ios::left);
        oss.str("");
        oss<<res->getName()<<res->getFNum();
        cout<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) cout<<" ";
        cout.width(4); cout.unsetf(ios::left);cout.setf(ios::right);
        if (res->getRootAtom() == (Atom*)NULL || mole->getMoleType()!=MoleType::PROTEIN)
        {
            // No Root Atom => So we put the first one in the list
            for (ItCAtom itA = res->firstAtom();itA != res->lastAtom();++itA)
            {
                if (AtmAccess[(*itA)->getNum()] != 0)
                {cout<< AtmAccess[(*itA)->getNum()]<<" ";break;}
            }

        }
        else cout<< AtmAccess[res->getRootAtom()->getNum()]<<" ";

        switch (res->getResType())
        {
            case ResType::MOD_AA:
            case ResType::STD_AA:
            case ResType::WATER:cout<<"RESIDUE ";break;
            default :             cout<<"GROUP   ";break;
        }

        if (mole->getMoleType() != MoleType::PROTEIN
            && mole->getMoleType() != MoleType::INT_CENT
            && mole->getMoleType() != MoleType::INT_PROT
            && mole->getMoleType() != MoleType::INT_LIG
            && mole->getMoleType() != MoleType::INT_MERG){cout<<"\n";continue;}
        cout << "     1 ";
        cout.width(4);
        cout.unsetf(ios::right);cout.setf(ios::left);
        if (res->getChain() != (Chain*)NULL) cout << res->getChainName()<<" ";else cout<< "X ";
        cout.width(4);
        cout <<res->getName()<<" ";
        cout.width(5); cout<< res->getNInterResidu()<<" ";
        if (currChain != res->getChainName() && res->getNInterResidu() != 2){cout <<"ROOT";currChain=res->getChainName();}
        else cout<< "****";
        cout<<" "<< res->getName()<<" "<<currChain<<" ";


        cout << res->getNum();
        cout<<"\n";

    }
    if (maxUnkAtom)
    {
        AccVal=ResAccess[mole->getUnknownRes().getNum()];
        cout.setf(ios::right);  cout.width(6);       cout<<AccVal<<" ";
        cout.unsetf(ios::right);cout.setf(ios::left);
        oss.str("");
        oss<<"UNKR"<<mole->getUnknownRes().getFNum();
        cout<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) cout<<" ";
        cout.width(4); cout.unsetf(ios::left);cout.setf(ios::right);
        for (ItCAtom itA = mole->getUnknownRes().firstAtom();itA != mole->getUnknownRes().lastAtom();++itA)
        {
            if (AtmAccess[(*itA)->getNum()] != 0)
            {cout<< AtmAccess[(*itA)->getNum()]<<" ";break;}
        }
        cout<<"GROUP   ";
        cout<<"\n";




    }


    if (!withSets && !withCycles&& NSet!=0) {cout<<endl; return;}
    cout<< "@<TRIPOS>SET"<<endl;
    unsigned int SetId=0;
    if (withSets)
    {

        // Starting with Atom Sets :
        for(std::vector<AtomSet>::const_iterator itAS = mole->firstAtomSet(); itAS != mole->lastAtomSet(); itAS++)
        {
            cout<< (*itAS).name<<" STATIC     ATOMS     "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique(const_cast<AtomList&>((*itAS).list));
            for (ItCAtom itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                const Atom& atm=**itAA;
                cout <<fname<<"\tATM\t"<< atm.getIdentifier()<<" " << atm.getNum()<< " " << AtmAccess[atm.getNum()]<<endl;
                if (AtmAccess[atm.getNum()]==0)continue;
                oss<< AtmAccess[atm.getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            cout << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }

        for(std::vector<ResSet>::iterator itAS = mole->firstResSet(); itAS != mole->lastResSet(); itAS++)
        {
            cout<< (*itAS).name<<" STATIC     SUBSTS    "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique((*itAS).list);
            for (ItRes itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                Residu &res = **itAA;
                if (ResAccess[res.getNum()]==0)continue;
                if (AtmAccess[res.getRootAtom()->getNum()] == 0)continue;
                cout << fname<<"\tRES\t"<<res.getRootAtom()->getIdentifier()<<" " << res.getRootAtom()->getNum()<< " " << AtmAccess[res.getRootAtom()->getNum()]<<endl;
                oss<< AtmAccess[res.getRootAtom()->getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            cout << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }


        if (ossCol.str().length()) cout <<ossCol.str();

    }
    cout<<endl;
}

/**
 * @brief MoleWriter::writeMOL2
 * @param mole : molecule to save
 *
 * save the given molecule into the specified file
 */
void MoleWriter::writeMOL2( Molecule * const mole, const bool outints)
{
    if (!ofstr.is_open()) {cerr << "File not opened" <<endl;return;}
    if (mole== (Molecule*)NULL)return;

    time_t timestamp = time(NULL);
    char buffer[256];
    unsigned int Natm=0;
    unsigned int  Nbd=0;
    unsigned int NRes=0;
    unsigned int NSet=0;
    unsigned int  RId=0;
    int RNum=0;
    unsigned int AccVal=0;
    string RName="";

    const size_t maxAtm= mole->getMaxNumAtom();
    const size_t maxCyclAtom=mole->numCyclAtom();
    const size_t maxAromAtom=mole->numAromAtom();
    const size_t maxUnkAtom =mole->numUnkAtom();
    const size_t nAllAtom = maxAtm+maxCyclAtom+maxAromAtom+maxUnkAtom;

    size_t AtmAccess[nAllAtom];
    size_t BdAccess[mole->getMaxNumBond()];
    size_t ResAccess[mole->getMaxNumResidu()+3];
    Coords coo;


    mole->getCycleRes().setNum(mole->getMaxNumResidu()+1);
    mole->getAromaticRes().setNum(mole->getMaxNumResidu()+2);
    mole->getUnknownRes().setNum(mole->getMaxNumResidu()+3);


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// INITIALIZATION ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    for (size_t i=0;i<nAllAtom;i++){AtmAccess[i]=0;   }
    for (size_t i=0;i<mole->getMaxNumBond();i++){BdAccess[i]=0;}
    for (size_t i=0;i<mole->numResidus()+3;i++){ResAccess[i]=0;}


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////// PRE-CALCULATIONS /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;

        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {

            Atom *atmNP = *itA;

            if (!withH && atmNP->isHydrogen()){AtmAccess[atmNP->getNum()]=0;continue; }
            if (onlyUsed && !atmNP->isUsed()){AtmAccess[atmNP->getNum()]=0;continue; }
            Natm++;
            //cout << atmNP->getIdentifier()<<" " <<atmNP->getNum()<<" " << Natm<<endl;
            AtmAccess[atmNP->getNum()]=Natm;
        }
    }

    for (ItCBond itB = mole->firstBond();
         itB != mole->lastBond();
         itB++)
    {
        try{
            if ((AtmAccess[(*itB)->getAtom1().getNum()] == 0
                 ||AtmAccess[(*itB)->getAtom2().getNum()] == 0)
                    || (onlyUsed && !(*itB)->isUsed()))
            {
                BdAccess[(*itB)->getNum()]=0;
                continue;
            }
        }
        catch (...)
        {
            cerr << "ERRS : "<< (*itB)->getAtom1().getIdentifier()<<"\t"<<(*itB)->getAtom2().getIdentifier()<<endl;
        }
        Nbd++;
        BdAccess[(*itB)->getNum()]=Nbd;
    }
    bool auth=false;
    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {
        auth=false;
        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {

            if (AtmAccess[(*itA)->getNum()] != 0){auth=true;break;}
        }
        if (!auth) {ResAccess[res->getNum()]=0;continue;}
        NRes++;
        ResAccess[res->getNum()]=NRes;


    }
    if (maxUnkAtom){NRes++;ResAccess[mole->getCycleRes().getNum()]=NRes;}
    if (withCycles)
    {
        if (maxCyclAtom >0){NRes++;ResAccess[mole->getCycleRes().getNum()]=NRes;}
        if (maxAromAtom >0){NRes++;ResAccess[mole->getAromaticRes().getNum()]=NRes;}
        NSet+=maxCyclAtom+maxAromAtom;
    }

    ostringstream ossCol;
    if (withSets)  { NSet+=mole->getNumSet();

        vector<unsigned int> TablCol[NB_COLOR];
        for(ItCAtom itAS = mole->firstColAtom(); itAS != mole->lastColAtom(); itAS++)
        {
            const Atom& atm = **itAS;
            if (AtmAccess[atm.getNum()]==0)continue;


            for (unsigned int ICol = 0; ICol < NB_COLOR;ICol++)
            {

                if (atm.getRcolor() == ColorData[ICol].RED
                        && atm.getGcolor() == ColorData[ICol].GREEN
                        && atm.getBcolor() == ColorData[ICol].BLUE)
                {

                    TablCol[ICol].push_back(AtmAccess[atm.getNum()]);break;}
            }
        }

        for (unsigned int ICol=0; ICol < NB_COLOR; ICol++)
        {

            if (TablCol[ICol].size()== 0) continue;
            NSet++;
            ossCol<< "ATOM$"<<ColorData[ICol].name<<" STATIC     ATOMS   COLORGROUP SYSTEM "<<endl;
            ossCol<< TablCol[ICol].size()<<" ";
            unsigned int count=1,countTOT=0;
            for (vector<unsigned int>::iterator it = TablCol[ICol].begin(); it != TablCol[ICol].end(); it++)
            {
                ossCol<< *it<<" ";count++;countTOT++;
                if (count == 18&& countTOT != TablCol[ICol].size()) {count=0; ossCol <<"\\"<<endl;}
            }
            ossCol<<endl;
        }


    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// HEADER BLOCK ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////


    ofstr<< "#    Name: "<< mole->getName()<<"\n";
    strftime(buffer, sizeof(buffer), "%c", localtime(&timestamp));
    ofstr << "#    Modification time : " << buffer<<"\n"<<"\n";
    ofstr <<"@<TRIPOS>MOLECULE"<<"\n";
    ofstr<< mole->getName()<<"\n";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<Natm<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<Nbd<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<NRes<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<"0"<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<NSet<<"\n";

    // Molecular name :
    switch (mole->getMoleType())
    {
    case MoleType::PROTEIN: ofstr<< "PROTEIN"<<"\n";break;
    case MoleType::COFACTOR:
    case MoleType::LIGAND:
    case MoleType::WATER:
    default: ofstr<<"SMALL"<<"\n";break;
    }

    ofstr << "NO_CHARGES"<<"\n"<<"\n"<<"\n";
    if (mole->getMoleType() == MoleType::PROTEIN)
        ofstr<< "@<TRIPOS>DICT"<<"\n"<<"PROTEIN macromol"<<"\n";


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// ATOM BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    ofstr<< "@<TRIPOS>ATOM"<<"\n";
    ofstr.setf(ios_base::fixed,ios_base::floatfield);ofstr.precision(4);
    unsigned int IdAtom,IdBond;
    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {
            Atom *atm=*itA;
            IdAtom=AtmAccess[atm->getNum()];
            if (IdAtom==0 )continue;

            RName=atm->getResidu()->getName();
            RId=ResAccess[atm->getResidu()->getNum()];
            RNum=atm->getResidu()->getFNum();


            if (withRotPos) coo = atm->rotpos; else coo = atm->fixpos;

            ofstr.setf(ios::right);ofstr.width(7);ofstr<<IdAtom<<" ";
            ofstr.unsetf(ios::right);ofstr.setf(ios::left); ofstr.width(7);ofstr<<atm->getName()<<" ";
            ofstr.unsetf(ios::left); ofstr.setf(ios::right);
            ofstr.width(10);ofstr<<coo.x<<" ";
            ofstr.width(9); ofstr<<coo.y<<" ";
            ofstr.width(9); ofstr<<coo.z<<" ";
            ofstr.width(6); ofstr.unsetf(ios::right);ofstr.setf(ios::left);  ofstr<<atm->getMOL2Type()<<" ";
            ofstr.width(4); ofstr.unsetf(ios::left); ofstr.setf(ios::right); ofstr<<RId<<" ";
            if (RName.length() >=3)
            {
                ofstr.width(3); ofstr.unsetf(ios::right);ofstr.setf(ios::left);   ofstr<<RName;
                ofstr.width(5);                                                  ofstr<<RNum<<" ";
            }
            else if (RName.length()==2)
            {
                ofstr.width(2); ofstr.unsetf(ios::right);ofstr.setf(ios::left);   ofstr<<RName;
                ofstr.width(6);                                                  ofstr<<RNum<<" ";

            }
            ofstr.width(9); ofstr.unsetf(ios::left); ofstr.setf(ios::right); ofstr<<atm->getPartialCharge()<<" ";
            if (atm->getBits().length()) ofstr << atm->getBits();
            ofstr<<"\n";
        }}

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    ofstr<< "@<TRIPOS>BOND"<<"\n";
    for (ItCBond it = mole->firstBond(); it != mole->lastBond();it++)
    {
        Bond *bd=*it;
        IdBond=BdAccess[bd->getNum()];
        if (IdBond == 0) continue;
        ofstr.setf(ios::right);ofstr.width(6);ofstr<<IdBond<<" ";
        ofstr.width(4); ofstr<<AtmAccess[bd->getAtom1().getNum()]<<" ";
        ofstr.width(4); ofstr<<AtmAccess[bd->getAtom2().getNum()]<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);ofstr.width(4);
        switch (bd->getBondType())
        {
        case BondType::SINGLE: ofstr<<"1";break;
        case BondType::DOUBLE: ofstr<<"2";break;
        case BondType::TRIPLE: ofstr<<"3";break;
        case BondType::AMIDE:  ofstr<<"am";break;
        case BondType::AROMATIC:  ofstr<<"ar";break;
        case BondType::ANY:  ofstr<<"an";break;
        case BondType::DUMMY:  ofstr<<"du";break;
        default:  ofstr<<"un";break;
        }
        ofstr<<" ";
        if (bd->getBits().length()) ofstr << bd->getBits();
        ofstr<<"\n";
    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// INTS BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    if (outints){
        ofstr<< "@<TRIPOS>INTS"<<"\n";
        ofstr.setf(ios_base::fixed,ios_base::floatfield);ofstr.precision(4);


        for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
        {

            Residu *res = *it;
            for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
            {
                Atom *atm=*itA;
                IdAtom=AtmAccess[atm->getNum()];
                if (IdAtom==0 )continue;
                if (withRotPos) coo = atm->rotpos; else coo = atm->fixpos;

                ofstr.setf(ios::right);ofstr.width(7);ofstr<<IdAtom<<" ";
                ofstr.setf(ios::right);ofstr.width(7);ofstr<<atm->getTNum();
//                for (int ints_ref : atm->getTNum()){
//                    ofstr.setf(ios::right);ofstr.width(7);ofstr<<ints_ref <<" ";
//
//                }
                ofstr <<"\n";
            }

        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////// SUBSTRUCTURE BLOCK ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    string currChain="";

    ofstr<< "@<TRIPOS>SUBSTRUCTURE"<<"\n";
    ostringstream oss;
    for (ItCRes it = mole->firstResidu(); it != mole->lastResidu();it++)
    {
        Residu *res = *it;
        AccVal=ResAccess[res->getNum()];
        if (AccVal == 0) continue;
        ofstr.setf(ios::right);ofstr.width(6);ofstr<<AccVal<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        oss.str("");
        oss<<res->getName()<<res->getFNum();
        ofstr<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) ofstr<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left);ofstr.setf(ios::right);
        if (res->getRootAtom() == (Atom*)NULL || mole->getMoleType()!=MoleType::PROTEIN)
        {
            // No Root Atom => So we put the first one in the list
            for (ItCAtom itA = res->firstAtom();itA != res->lastAtom();++itA)
            {
                if (AtmAccess[(*itA)->getNum()] != 0)
                {ofstr<< AtmAccess[(*itA)->getNum()]<<" ";break;}
            }

        }
        else ofstr<< AtmAccess[res->getRootAtom()->getNum()]<<" ";

        switch (res->getResType())
        {
        case ResType::MOD_AA:
        case ResType::STD_AA:
        case ResType::WATER:ofstr<<"RESIDUE ";break;
        default :             ofstr<<"GROUP   ";break;
        }

        if (mole->getMoleType() != MoleType::PROTEIN
                && mole->getMoleType() != MoleType::INT_CENT
                && mole->getMoleType() != MoleType::INT_PROT
                && mole->getMoleType() != MoleType::INT_LIG
                && mole->getMoleType() != MoleType::INT_MERG){ofstr<<"\n";continue;}
        ofstr << "     1 ";
        ofstr.width(4);
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        if (res->getChain() != (Chain*)NULL) ofstr << res->getChainName()<<" ";else ofstr<< "X ";
        ofstr.width(4);
        ofstr <<res->getName()<<" ";
        ofstr.width(5); ofstr<< res->getNInterResidu()<<" ";
        if (currChain != res->getChainName() && res->getNInterResidu() != 2){ofstr <<"ROOT";currChain=res->getChainName();}
        else ofstr<< "****";
        ofstr<<" "<< res->getName()<<" "<<currChain<<" ";


        ofstr << res->getNum();
        ofstr<<"\n";

    }
    if (maxUnkAtom)
    {
        cout << "ATOM IN MAXUNK<<"<<"\n";
        AccVal=ResAccess[mole->getUnknownRes().getNum()];
        ofstr.setf(ios::right);  ofstr.width(6);       ofstr<<AccVal<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        oss.str("");
        oss<<"UNK"<<mole->getUnknownRes().getFNum();
        ofstr<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) ofstr<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left);ofstr.setf(ios::right);
        for (ItCAtom itA = mole->getUnknownRes().firstAtom();itA != mole->getUnknownRes().lastAtom();++itA)
        {
            if (AtmAccess[(*itA)->getNum()] != 0)
            {ofstr<< AtmAccess[(*itA)->getNum()]<<" ";break;}
        }
        ofstr<<"GROUP   ";
        ofstr<<"\n";




    }


    if (!withSets && !withCycles&& NSet!=0) {ofstr<<endl; return;}
    ofstr<< "@<TRIPOS>SET"<<endl;
    unsigned int SetId=0;
    if (withSets)
    {

        // Starting with Atom Sets :
        for(std::vector<AtomSet>::const_iterator itAS = mole->firstAtomSet(); itAS != mole->lastAtomSet(); itAS++)
        {
            ofstr<< (*itAS).name<<" STATIC     ATOMS     "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique(const_cast<AtomList&>((*itAS).list));
            for (ItCAtom itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                const Atom& atm=**itAA;
                cout <<fname<<"\tATM\t"<< atm.getIdentifier()<<" " << atm.getNum()<< " " << AtmAccess[atm.getNum()]<<endl;
                if (AtmAccess[atm.getNum()]==0)continue;
                oss<< AtmAccess[atm.getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            ofstr << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }

        for(std::vector<ResSet>::iterator itAS = mole->firstResSet(); itAS != mole->lastResSet(); itAS++)
        {
            ofstr<< (*itAS).name<<" STATIC     SUBSTS    "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique((*itAS).list);
            for (ItRes itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                Residu &res = **itAA;
                if (ResAccess[res.getNum()]==0)continue;
                if (AtmAccess[res.getRootAtom()->getNum()] == 0)continue;
                cout << fname<<"\tRES\t"<<res.getRootAtom()->getIdentifier()<<" " << res.getRootAtom()->getNum()<< " " << AtmAccess[res.getRootAtom()->getNum()]<<endl;
                oss<< AtmAccess[res.getRootAtom()->getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            ofstr << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }


        if (ossCol.str().length()) ofstr <<ossCol.str();

    }

    /* if (withCycles)
    {
      unsigned int ncy=0;
      ostringstream ossCentroid;
      for (ICPIter<Cycle,Molecule> cy(mole);cy;++cy)
      {
        ncy++;
        SetId++;
        oss.str("");
        oss << "CYCLE_"<< ncy<< ((cy->isAromatic())?"_AROM":"")<<" STATIC     ATOMS     "<<endl ;
        bool ok=true;
        oss<< cy->getNumAtom()<<" ";
        for (ItCAtom itA = cy->first(); itA != cy->end(); itA++)
          {
            ITTATM = AtmAccess.find(const_cast<Atom*>(*itA));
            if (ITTATM == AtmAccess.end()){ok=false;break;}
            oss << (*ITTATM).second<<" ";
          }
        oss<<endl;

        if (ok) ofstr << oss.str();
        ossCentroid<< "CENTR_CYCLE"<<ncy<<endl<<AtmAccess.at(cy->getCenter())<<" "<< SetId<<endl;
      }
      if (ossCentroid.str().length()) ofstr << "@<TRIPOS>CENTROID"<<endl<< ossCentroid.str()<<endl;


    }*/

    ofstr<<endl;
}


/**
 * @brief MoleWriter::writePharmaLig
 * @param mole : molecule to save
 *
 * save the given molecule with pharmacophoric property into the specified file
 */
void MoleWriter::writePharmaLig( Molecule * const mole)
{
    if (!ofstr.is_open()) {cerr << "File not opened" <<endl;return;}
    if (mole== (Molecule*)NULL)return;

    time_t timestamp = time(NULL);
    char buffer[256];
    unsigned int Natm=0;
    unsigned int  Nbd=0;
    unsigned int NRes=0;
    unsigned int NSet=0;
    unsigned int  RId=0;
    int RNum=0;
    unsigned int AccVal=0;
    string RName="";

    const size_t maxAtm= mole->getMaxNumAtom();
    const size_t maxCyclAtom=mole->numCyclAtom();
    const size_t maxAromAtom=mole->numAromAtom();
    const size_t maxUnkAtom =mole->numUnkAtom();
    const size_t nAllAtom = maxAtm+maxCyclAtom+maxAromAtom+maxUnkAtom;

    size_t AtmAccess[nAllAtom];
    size_t BdAccess[mole->getMaxNumBond()];
    size_t ResAccess[mole->getMaxNumResidu()+3];
    Coords coo;


    mole->getCycleRes().setNum(mole->getMaxNumResidu()+1);
    mole->getAromaticRes().setNum(mole->getMaxNumResidu()+2);
    mole->getUnknownRes().setNum(mole->getMaxNumResidu()+3);


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// INITIALIZATION ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    for (size_t i=0;i<nAllAtom;i++){AtmAccess[i]=0;   }
    for (size_t i=0;i<mole->getMaxNumBond();i++){BdAccess[i]=0;}
    for (size_t i=0;i<mole->numResidus()+3;i++){ResAccess[i]=0;}


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////// PRE-CALCULATIONS /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;

        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {

            Atom *atmNP = *itA;

            if (!withH && atmNP->isHydrogen()){AtmAccess[atmNP->getNum()]=0;continue; }
            if (onlyUsed && !atmNP->isUsed()){AtmAccess[atmNP->getNum()]=0;continue; }
            Natm++;
            //cout << atmNP->getIdentifier()<<" " <<atmNP->getNum()<<" " << Natm<<endl;
            AtmAccess[atmNP->getNum()]=Natm;
        }
    }

    for (ItCBond itB = mole->firstBond();
         itB != mole->lastBond();
         itB++)
    {
        try{
            if ((AtmAccess[(*itB)->getAtom1().getNum()] == 0
                 ||AtmAccess[(*itB)->getAtom2().getNum()] == 0)
                    || (onlyUsed && !(*itB)->isUsed()))
            {
                BdAccess[(*itB)->getNum()]=0;
                continue;
            }
        }
        catch (...)
        {
            cerr << "ERRS : "<< (*itB)->getAtom1().getIdentifier()<<"\t"<<(*itB)->getAtom2().getIdentifier()<<endl;
        }
        Nbd++;
        BdAccess[(*itB)->getNum()]=Nbd;
    }
    bool auth=false;
    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {
        auth=false;
        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {

            if (AtmAccess[(*itA)->getNum()] != 0){auth=true;break;}
        }
        if (!auth) {ResAccess[res->getNum()]=0;continue;}
        NRes++;
        ResAccess[res->getNum()]=NRes;


    }
    if (maxUnkAtom){NRes++;ResAccess[mole->getCycleRes().getNum()]=NRes;}
    if (withCycles)
    {
        if (maxCyclAtom >0){NRes++;ResAccess[mole->getCycleRes().getNum()]=NRes;}
        if (maxAromAtom >0){NRes++;ResAccess[mole->getAromaticRes().getNum()]=NRes;}
        NSet+=maxCyclAtom+maxAromAtom;
    }

    ostringstream ossCol;
    if (withSets)  { NSet+=mole->getNumSet();

        vector<unsigned int> TablCol[NB_COLOR];
        for(ItCAtom itAS = mole->firstColAtom(); itAS != mole->lastColAtom(); itAS++)
        {
            const Atom& atm = **itAS;
            if (AtmAccess[atm.getNum()]==0)continue;


            for (unsigned int ICol = 0; ICol < NB_COLOR;ICol++)
            {

                if (atm.getRcolor() == ColorData[ICol].RED
                        && atm.getGcolor() == ColorData[ICol].GREEN
                        && atm.getBcolor() == ColorData[ICol].BLUE)
                {

                    TablCol[ICol].push_back(AtmAccess[atm.getNum()]);break;}
            }
        }

        for (unsigned int ICol=0; ICol < NB_COLOR; ICol++)
        {

            if (TablCol[ICol].size()== 0) continue;
            NSet++;
            ossCol<< "ATOM$"<<ColorData[ICol].name<<" STATIC     ATOMS   COLORGROUP SYSTEM "<<endl;
            ossCol<< TablCol[ICol].size()<<" ";
            unsigned int count=1,countTOT=0;
            for (vector<unsigned int>::iterator it = TablCol[ICol].begin(); it != TablCol[ICol].end(); it++)
            {
                ossCol<< *it<<" ";count++;countTOT++;
                if (count == 18&& countTOT != TablCol[ICol].size()) {count=0; ossCol <<"\\"<<endl;}
            }
            ossCol<<endl;
        }


    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// HEADER BLOCK ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////


    ofstr<< "#    Name: "<< mole->getName()<<"\n";
    strftime(buffer, sizeof(buffer), "%c", localtime(&timestamp));
    ofstr << "#    Modification time : " << buffer<<"\n"<<"\n";
    ofstr <<"@<TRIPOS>MOLECULE"<<"\n";
    ofstr<< mole->getName()<<"\n";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<Natm<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<Nbd<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<NRes<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<"0"<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<NSet<<"\n";

    // Molecular name :
    switch (mole->getMoleType())
    {
    case MoleType::PROTEIN: ofstr<< "PROTEIN"<<"\n";break;
    case MoleType::COFACTOR:
    case MoleType::LIGAND:
    case MoleType::WATER:
    default: ofstr<<"SMALL"<<"\n";break;
    }

    ofstr << "NO_CHARGES"<<"\n"<<"\n"<<"\n";
    if (mole->getMoleType() == MoleType::PROTEIN)
        ofstr<< "@<TRIPOS>DICT"<<"\n"<<"BIOPOLYMER macromol"<<"\n";


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// ATOM BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    ofstr<< "@<TRIPOS>ATOM"<<"\n";
    ofstr.setf(ios_base::fixed,ios_base::floatfield);ofstr.precision(4);
    unsigned int IdAtom,IdBond;
    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {
            Atom *atm=*itA;
            IdAtom=AtmAccess[atm->getNum()];
            if (IdAtom==0 )continue;
            string namepharma = "";
            if (atm->props.isAcceptor() && atm->props.isDonor() ){
                namepharma="OG";
            }
            else if (atm->props.isAcceptor()){namepharma="O";}
            else if (atm->props.isDonor()){namepharma="N";}
            else if (atm->props.isAnion()){namepharma="OD1";}
            else if (atm->props.isCation()){namepharma="NZ";}
            else if (atm->props.isHydrophobic()){namepharma="CA";}
            else if (atm->props.isAromatic()){namepharma="CZ";}
            else { namepharma=atm->getName();}

            RName=atm->getResidu()->getName();
            RId=ResAccess[atm->getResidu()->getNum()];
            RNum=atm->getResidu()->getFNum();


            if (withRotPos) coo = atm->rotpos; else coo = atm->fixpos;

            ofstr.setf(ios::right);ofstr.width(7);ofstr<<IdAtom<<" ";
            ofstr.unsetf(ios::right);ofstr.setf(ios::left); ofstr.width(7);ofstr<<namepharma<<" ";
            ofstr.unsetf(ios::left); ofstr.setf(ios::right);
            ofstr.width(10);ofstr<<coo.x<<" ";
            ofstr.width(9); ofstr<<coo.y<<" ";
            ofstr.width(9); ofstr<<coo.z<<" ";
            ofstr.width(6); ofstr.unsetf(ios::right);ofstr.setf(ios::left);  ofstr<<atm->getMOL2Type()<<" ";
            ofstr.width(4); ofstr.unsetf(ios::left); ofstr.setf(ios::right); ofstr<<RId<<" ";
            if (RName.length() >=3)
            {
                ofstr.width(3); ofstr.unsetf(ios::right);ofstr.setf(ios::left);   ofstr<<RName;
                ofstr.width(5);                                                  ofstr<<RNum<<" ";
            }
            else if (RName.length()==2)
            {
                ofstr.width(2); ofstr.unsetf(ios::right);ofstr.setf(ios::left);   ofstr<<RName;
                ofstr.width(6);                                                  ofstr<<RNum<<" ";

            }
            ofstr.width(9); ofstr.unsetf(ios::left); ofstr.setf(ios::right); ofstr<<atm->getPartialCharge()<<" ";
            if (atm->getBits().length()) ofstr << atm->getBits();
            ofstr<<"\n";
        }}

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    ofstr<< "@<TRIPOS>BOND"<<"\n";
    for (ItCBond it = mole->firstBond(); it != mole->lastBond();it++)
    {
        Bond *bd=*it;
        IdBond=BdAccess[bd->getNum()];
        if (IdBond == 0) continue;
        ofstr.setf(ios::right);ofstr.width(6);ofstr<<IdBond<<" ";
        ofstr.width(4); ofstr<<AtmAccess[bd->getAtom1().getNum()]<<" ";
        ofstr.width(4); ofstr<<AtmAccess[bd->getAtom2().getNum()]<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);ofstr.width(4);
        switch (bd->getBondType())
        {
        case BondType::SINGLE: ofstr<<"1";break;
        case BondType::DOUBLE: ofstr<<"2";break;
        case BondType::TRIPLE: ofstr<<"3";break;
        case BondType::AMIDE:  ofstr<<"am";break;
        case BondType::AROMATIC:  ofstr<<"ar";break;
        case BondType::ANY:  ofstr<<"an";break;
        case BondType::DUMMY:  ofstr<<"du";break;
        default:  ofstr<<"un";break;
        }
        ofstr<<" ";
        if (bd->getBits().length()) ofstr << bd->getBits();
        ofstr<<"\n";
    }


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    string currChain="";

    ofstr<< "@<TRIPOS>SUBSTRUCTURE"<<"\n";
    ostringstream oss;
    for (ItCRes it = mole->firstResidu(); it != mole->lastResidu();it++)
    {
        Residu *res = *it;
        AccVal=ResAccess[res->getNum()];
        if (AccVal == 0) continue;
        ofstr.setf(ios::right);ofstr.width(6);ofstr<<AccVal<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        oss.str("");
        oss<<res->getName()<<res->getFNum();
        ofstr<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) ofstr<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left);ofstr.setf(ios::right);
        if (res->getRootAtom() == (Atom*)NULL || mole->getMoleType()!=MoleType::PROTEIN)
        {
            // No Root Atom => So we put the first one in the list
            for (ItCAtom itA = res->firstAtom();itA != res->lastAtom();++itA)
            {
                if (AtmAccess[(*itA)->getNum()] != 0)
                {ofstr<< AtmAccess[(*itA)->getNum()]<<" ";break;}
            }

        }
        else ofstr<< AtmAccess[res->getRootAtom()->getNum()]<<" ";

        switch (res->getResType())
        {
        case ResType::MOD_AA:
        case ResType::STD_AA:
        case ResType::WATER:ofstr<<"RESIDUE ";break;
        default :             ofstr<<"GROUP   ";break;
        }

        if (mole->getMoleType() != MoleType::PROTEIN
                && mole->getMoleType() != MoleType::INT_CENT
                && mole->getMoleType() != MoleType::INT_PROT
                && mole->getMoleType() != MoleType::INT_LIG
                && mole->getMoleType() != MoleType::INT_MERG){ofstr<<"\n";continue;}
        ofstr << "     1 ";
        ofstr.width(4);
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        if (res->getChain() != (Chain*)NULL) ofstr << res->getChainName()<<" ";else ofstr<< "X ";
        ofstr.width(4);
        ofstr <<res->getName()<<" ";
        ofstr.width(5); ofstr<< res->getNInterResidu()<<" ";
        if (currChain != res->getChainName() && res->getNInterResidu() != 2){ofstr <<"ROOT";currChain=res->getChainName();}
        else ofstr<< "****";
        ofstr<<" "<< res->getName()<<" "<<currChain<<" ";


        ofstr << res->getNum();
        ofstr<<"\n";

    }
    if (maxUnkAtom)
    {
        cout << "ATOM IN MAXUNK<<"<<"\n";
        AccVal=ResAccess[mole->getUnknownRes().getNum()];
        ofstr.setf(ios::right);  ofstr.width(6);       ofstr<<AccVal<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        oss.str("");
        oss<<"UNK"<<mole->getUnknownRes().getFNum();
        ofstr<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) ofstr<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left);ofstr.setf(ios::right);
        for (ItCAtom itA = mole->getUnknownRes().firstAtom();itA != mole->getUnknownRes().lastAtom();++itA)
        {
            if (AtmAccess[(*itA)->getNum()] != 0)
            {ofstr<< AtmAccess[(*itA)->getNum()]<<" ";break;}
        }
        ofstr<<"GROUP   ";
        ofstr<<"\n";




    }


    if (!withSets && !withCycles&& NSet!=0) {ofstr<<endl; return;}
    ofstr<< "@<TRIPOS>SET"<<endl;
    unsigned int SetId=0;
    if (withSets)
    {

        // Starting with Atom Sets :
        for(std::vector<AtomSet>::const_iterator itAS = mole->firstAtomSet(); itAS != mole->lastAtomSet(); itAS++)
        {
            ofstr<< (*itAS).name<<" STATIC     ATOMS     "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique(const_cast<AtomList&>((*itAS).list));
            for (ItCAtom itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                const Atom& atm=**itAA;
                cout <<fname<<"\tATM\t"<< atm.getIdentifier()<<" " << atm.getNum()<< " " << AtmAccess[atm.getNum()]<<endl;
                if (AtmAccess[atm.getNum()]==0)continue;
                oss<< AtmAccess[atm.getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            ofstr << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }

        for(std::vector<ResSet>::iterator itAS = mole->firstResSet(); itAS != mole->lastResSet(); itAS++)
        {
            ofstr<< (*itAS).name<<" STATIC     SUBSTS    "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique((*itAS).list);
            for (ItRes itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                Residu &res = **itAA;
                if (ResAccess[res.getNum()]==0)continue;
                if (AtmAccess[res.getRootAtom()->getNum()] == 0)continue;
                cout << fname<<"\tRES\t"<<res.getRootAtom()->getIdentifier()<<" " << res.getRootAtom()->getNum()<< " " << AtmAccess[res.getRootAtom()->getNum()]<<endl;
                oss<< AtmAccess[res.getRootAtom()->getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            ofstr << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }


        if (ossCol.str().length()) ofstr <<ossCol.str();

    }

    /* if (withCycles)
    {
      unsigned int ncy=0;
      ostringstream ossCentroid;
      for (ICPIter<Cycle,Molecule> cy(mole);cy;++cy)
      {
        ncy++;
        SetId++;
        oss.str("");
        oss << "CYCLE_"<< ncy<< ((cy->isAromatic())?"_AROM":"")<<" STATIC     ATOMS     "<<endl ;
        bool ok=true;
        oss<< cy->getNumAtom()<<" ";
        for (ItCAtom itA = cy->first(); itA != cy->end(); itA++)
          {
            ITTATM = AtmAccess.find(const_cast<Atom*>(*itA));
            if (ITTATM == AtmAccess.end()){ok=false;break;}
            oss << (*ITTATM).second<<" ";
          }
        oss<<endl;

        if (ok) ofstr << oss.str();
        ossCentroid<< "CENTR_CYCLE"<<ncy<<endl<<AtmAccess.at(cy->getCenter())<<" "<< SetId<<endl;
      }
      if (ossCentroid.str().length()) ofstr << "@<TRIPOS>CENTROID"<<endl<< ossCentroid.str()<<endl;


    }*/

    ofstr<<endl;
}



/**
 * @brief MoleWriter::writeMOL2
 * @param mole : molecule to save
 *
 * save substructure of the the given molecule into the specified file
 */
void MoleWriter::writeMOL2ss(Molecule * const mole , const int mod)
{
    double maxdist=0;
    double dist;
    // Feature importante, depart des pharmacophores
    Coords Eatm[6];
    Molecule Matm[6];
    string outname[6];

    for (int i =0;i<4;i++){
        Eatm[i]= Coords(0,0,0);
        ostringstream o;
        o << "Pharma_" << mod << "_ss" << i;
        Matm[i].setName(o.str());
        o << ".mol2";
        outname[i] = o.str();
    }

    int flag=0;
    double div=1;
    int nbfeat =0;

    for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
    {

        Residu *res = *it;
        for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
        {
            Atom *atm=*itA;

            if (atm->getMOL2Type() != "H" ){
                nbfeat++;
                continue;
            }

            for (ItCAtom itA2 = res->firstAtom(); itA2 != res->lastAtom(); itA2++)
            {
                Atom *atm2=*itA2;
//                cout << atm2->getAtomicName()<< endl;
//                if (atm2->getMOL2Type() != "H" ){ continue; }
                dist =atm->fixpos.calcDist(atm2->fixpos);
                if ( dist > maxdist){
                    maxdist = dist;
                    if (Eatm[0] != Coords(0,0,0)){
                        Coords temp1 = Eatm[0];
                        Coords temp2 = Eatm[1];
                        Eatm[3].setCoords(temp1);
                        Eatm[4].setCoords(temp2);
                    }
                    Eatm[0].setCoords(atm->fixpos);
                    Eatm[1].setCoords(atm2->fixpos);
                }
            }
        }
    }

//    cout << nbfeat << endl;
    if (nbfeat < 20){
        flag=0;
        div=1.0;
    }else if (nbfeat < 30){
        flag=2;
        div=0.6;

    }else if (nbfeat < 40){
        flag = 3;
        div = 0.5;
    }else {
        flag = 6;
        div = 0.5;
    }

    flag = 0;

    double div_temp = 0;
    Eatm[2].setCoords((Eatm[0].x+Eatm[1].x)/2,(Eatm[0].y+Eatm[1].y)/2,(Eatm[0].z+Eatm[1].z)/2);
    Eatm[5].setCoords((Eatm[3].x+Eatm[4].x)/2,(Eatm[3].y+Eatm[4].y)/2,(Eatm[3].z+Eatm[4].z)/2);
    int posi=0;

    for (posi=0 ; posi < flag; posi++){
        if (posi == 2 || posi == 5) {
            div_temp = div/2;
        }else {
            div_temp = div;
        }

        //    Molecule
        Residu &resSER = Matm[posi].addResidu("X",1,"SER");
        Residu &resALA = Matm[posi].addResidu("X",2,"ALA");
        Residu &resASP = Matm[posi].addResidu("X",3,"ASP");
        Residu &resLYS = Matm[posi].addResidu("X",4,"LYS");
        Residu &resPHE = Matm[posi].addResidu("X",5,"PHE");
        Residu &resGLY = Matm[posi].addResidu("X",6,"GLY");
        Residu &resTRP = Matm[posi].addResidu("X",7,"TRP");
        Residu &resCUB = Matm[posi].addResidu("X",8,"CUB");

        for (ItCRes it = mole->firstResidu();it != mole->lastResidu(); it++)
        {

            Residu *res = *it;
            for (ItCAtom itA = res->firstAtom(); itA != res->lastAtom(); itA++)
            {
                Atom *atm=*itA;
//                cout <<  Eatm[posi].calcDist(atm->fixpos) << " max " << maxdist*div_temp<<" ou " <<  atm->getMOL2Type()<<endl;
                if (Eatm[posi].calcDist(atm->fixpos) > (maxdist*div_temp) || atm->getMOL2Type() == "H"){
                    continue;}
//                cout << "ca passe " << atm->props.toString() << endl;
                if(atm->getName() == "Zn"){
                    Matm[posi].addAtom("Zn",atm->fixpos,"Zn","Zn", &resTRP);
                }
                if(atm->getName() == "OD1"){
                    Matm[posi].addAtom("O",atm->fixpos,"OD1","O.co2", &resASP);
                }
                if(atm->getName() == "NZ"){
                    Matm[posi].addAtom("N",atm->fixpos,"NZ","N.am" , &resLYS);
                }
                if(atm->getName() == "OG"){
                    Matm[posi].addAtom("O",atm->fixpos,"OG","O.3" , &resSER);
                }
                if(atm->getName() == "N"){
                    Matm[posi].addAtom("N",atm->fixpos,"N","N.am"  , &resALA);
                }
                if(atm->getName() == "O"){
                    Matm[posi].addAtom("O",atm->fixpos,"O","O.2"  , &resALA);
                }
                if(atm->getName() == "CA"){
                    Matm[posi].addAtom("C",atm->fixpos,"CA","C.3" , &resGLY);

                }
            }
        }
        //        cout << "Affichage protein" << Matm[posi].toString() << endl;
        cout << "Creation substructure : " << outname[posi] <<  endl;
        newFile(outname[posi],FileFormat::MOL2);
        writeMOL2(&Matm[posi]);
        // Mettre l'ecriture des fichier ici.

    }
}


void MoleWriter::writeMOL2( const Complex& complex, const std::string& name)
{
    if (!ofstr.is_open()) {cerr << "fILE NOT OPENED" <<endl;return;}

    vector<const Atom*> atomlist;
    vector<const Bond*> bondlist;
    vector<const Residu*> residulist;
    vector<const Chain*> chainlist;


    size_t mtype=0;
    for (mtype=0;mtype < NB_MOLETYPE; mtype++)
        for (ItCMole it = complex.firstMole(mtype);
             it !=complex.lastMole(mtype);
             it++)
        {
            const Molecule &mole=**it;
            for (ItCChain itCh=mole.firstChain();itCh!= mole.lastChain();++itCh)
            {
                const Chain& chain =**itCh;
                bool chainallowed=false;
                for (ItCRes itR =chain.firstC();itR != chain.lastC();++itR)
                {
                    const Residu& residu=**itR;
                    bool resAllowed=false;
                    for (ItCAtom itA = residu.firstCAtom();
                         itA!= residu.lastCAtom();
                         ++itA)
                    {
                        const Atom&atom =**itA;
                        if (onlyUsed && !atom.isUsed())continue;
                        if (!withH && atom.isHydrogen())continue;
                        resAllowed=true;
                        atomlist.push_back(&atom);
                    }
                    if (!resAllowed)continue;
                    residulist.push_back(&residu);
                    chainallowed=true;

                }
                if (!chainallowed)continue;
                chainlist.push_back(&chain);
            }
            for (ItCBond itB=mole.firstBond();
                 itB!=mole.lastBond();
                 ++itB)
            {
                const Bond &bond=**itB;
                if (onlyUsed && !bond.isUsed() ) continue;
                bondlist.push_back(&bond);
            }

        }

    unsigned int SetId=0;
    vector<AtomSet>  atomsets;
    if (withSets)
    { size_t mtype=0;
        for (;mtype < NB_MOLETYPE; mtype++)
            for (ItCMole it = complex.firstMole(mtype);
                 it !=complex.lastMole(mtype);
                 it++)
            {
                const Molecule &mole=**it;

                if (mole.getNumSet() >0)
                    for(std::vector<AtomSet>::const_iterator itAS = mole.firstAtomSet(); itAS != mole.lastAtomSet(); itAS++)
                    {
                        const AtomSet& moleset=*itAS;

                        bool found=false;
                        for (vector<AtomSet>::iterator itAsets=atomsets.begin(); itAsets !=atomsets.end();++itAsets)
                        {
                            if ((*itAsets).name==moleset.name)
                            {
                                found=true;
                                (*itAsets).list.insert((*itAsets).list.end(),moleset.list.begin(),moleset.list.end());
                            }
                        }
                        if (!found)
                        {
                            AtomSet newSet;
                            newSet.name=moleset.name;
                            newSet.subtype=moleset.subtype;
                            newSet.type=moleset.type;
                            newSet.list=moleset.list;
                            atomsets.push_back(newSet);
                        }
                    }
            }
    }


    const size_t maxAtm = atomlist.size();
    const size_t maxBond = bondlist.size();
    const size_t maxRes = residulist.size();
    const size_t maxChain= chainlist.size();
    time_t timestamp = time(NULL);
    char buffer[256];
    unsigned int Natm=0;
    unsigned int  Nbd=0;
    unsigned int NRes=0;
    unsigned int NSet=0;
    unsigned int  RId=0;



    ofstr<< "#    Name: "<< name<<"\n";
    strftime(buffer, sizeof(buffer), "%c", localtime(&timestamp));
    ofstr << "#    Modification time : " << buffer<<"\n"<<"\n";
    ofstr <<"@<TRIPOS>MOLECULE"<<"\n";
    ofstr<< name<<"\n";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<maxAtm<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<maxBond<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<maxRes<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<"0"<<" ";
    ofstr.setf(ios::right);ofstr.width(5); ofstr<<atomsets.size()<<"\n";
    ofstr<< "PROTEIN"<<"\n";
    ofstr << "NO_CHARGES"<<"\n"<<"\n"<<"\n";
    ofstr<< "@<TRIPOS>DICT"<<"\n"<<"BIOPOLYMER macromol"<<"\n";

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// ATOM BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    ofstr<< "@<TRIPOS>ATOM"<<"\n";
    ofstr.setf(ios_base::fixed,ios_base::floatfield);ofstr.precision(4);


    for (size_t iAtm=0; iAtm < maxAtm;++iAtm)
    {

        const Atom& atom = *atomlist.at(iAtm);

        const std::string& RName = atom.getResiduName();
        const size_t RId = std::distance(residulist.begin(),find(residulist.begin(),residulist.end(),atom.getResidu()))+1;
        const int RNum = atom.getResidu()->getFNum();

        ofstr.setf(ios::right);ofstr.width(7);ofstr<<(iAtm+1)<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left); ofstr.width(7);ofstr<<atom.getName()<<" ";
        ofstr.unsetf(ios::left); ofstr.setf(ios::right);
        if (withRotPos)
        {
            ofstr.width(10);ofstr<<atom.rotpos.x<<" ";
            ofstr.width(9); ofstr<<atom.rotpos.y<<" ";
            ofstr.width(9); ofstr<<atom.rotpos.z<<" ";
        }
        else
        {
            ofstr.width(10);ofstr<<atom.fixpos.x<<" ";
            ofstr.width(9); ofstr<<atom.fixpos.y<<" ";
            ofstr.width(9); ofstr<<atom.fixpos.z<<" ";
        }

        ofstr.width(6); ofstr.unsetf(ios::right);ofstr.setf(ios::left);  ofstr<<atom.getMOL2Type()<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left); ofstr.setf(ios::right); ofstr<<RId<<" ";
        if (RName.length() >=3)
        {
            ofstr.width(3); ofstr.unsetf(ios::right);ofstr.setf(ios::left);   ofstr<<RName;
            ofstr.width(5);                                                  ofstr<<RNum<<" ";
        }
        else if (RName.length()==2)
        {
            ofstr.width(2); ofstr.unsetf(ios::right);ofstr.setf(ios::left);   ofstr<<RName;
            ofstr.width(6);                                                  ofstr<<RNum<<" ";

        }
        ofstr.width(9); ofstr.unsetf(ios::left); ofstr.setf(ios::right); ofstr<<atom.getPartialCharge()<<" ";
        if (atom.getBits().length()) ofstr << atom.getBits();
        ofstr<<"\n";


    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    ofstr<< "@<TRIPOS>BOND"<<"\n";
    for (size_t iBd=0; iBd < maxBond; ++iBd)
    {
        const Bond &bd=*bondlist.at(iBd);
        ofstr.setf(ios::right);ofstr.width(6);ofstr<<(iBd+1)<<" ";
        ofstr.width(4); ofstr<< std::distance(atomlist.begin(),find(atomlist.begin(),atomlist.end(),&bd.getAtom1()))+1 <<" ";
        ofstr.width(4); ofstr<< std::distance(atomlist.begin(),find(atomlist.begin(),atomlist.end(),&bd.getAtom2()))+1<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);ofstr.width(4);
        switch (bd.getBondType())
        {
        case BondType::SINGLE: ofstr<<"1";break;
        case BondType::DOUBLE: ofstr<<"2";break;
        case BondType::TRIPLE: ofstr<<"3";break;
        case BondType::AMIDE:  ofstr<<"am";break;
        case BondType::AROMATIC:  ofstr<<"ar";break;
        case BondType::ANY:  ofstr<<"an";break;
        case BondType::DUMMY:  ofstr<<"du";break;
        default:  ofstr<<"un";break;
        }
        ofstr<<" ";
        if (bd.getBits().length()) ofstr << bd.getBits();
        ofstr<<"\n";
    }


    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// BOND BLOCK ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    string currChain="";

    ofstr<< "@<TRIPOS>SUBSTRUCTURE"<<"\n";
    ostringstream oss;
    for (size_t iRes=0; iRes < maxRes; ++iRes)
    {
        const Residu& residu = *residulist.at(iRes);

        ofstr.setf(ios::right);ofstr.width(6);ofstr<<(iRes+1)<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        oss.str("");
        oss<<residu.getName()<<residu.getFNum();
        ofstr<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) ofstr<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left);ofstr.setf(ios::right);
        if (residu.getRootAtom() == (Atom*)NULL)
        {
            // No Root Atom => So we put the first one in the list
            for (ItCAtom itA = residu.firstCAtom();itA != residu.lastCAtom();++itA)
            {
                const Atom& atm = **itA;
                const  vector<const Atom*>::iterator pos=find(atomlist.begin(),atomlist.end(),&atm);
                if (pos==atomlist.end())continue;
                ofstr<< std::distance(atomlist.begin(),pos)+1<<" ";break;
            }

        }
        else ofstr<< std::distance(atomlist.begin(),find(atomlist.begin(),atomlist.end(),residu.getRootAtom()))+1<<" ";

        switch (residu.getResType())
        {
        case ResType::MOD_AA:
        case ResType::STD_AA:
        case ResType::WATER:ofstr<<"RESIDUE ";break;
        default :             ofstr<<"GROUP   ";break;
        }

        ofstr << "     1 ";
        ofstr.width(4);
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        if (residu.getChain() != (Chain*)NULL) ofstr << residu.getChainName()<<" ";else ofstr<< "X ";
        ofstr.width(4);
        ofstr <<residu.getName()<<" ";
        ofstr.width(5); ofstr<< residu.getNInterResidu()<<" ";
        if (currChain != residu.getChainName() && residu.getNInterResidu() != 2){ofstr <<"ROOT";currChain=residu.getChainName();}
        else ofstr<< "****";
        ofstr<<" "<< residu.getName()<<" "<<currChain<<" ";


        ofstr << residu.getNum();
        ofstr<<"\n";

    }
    /* if (maxUnkAtom)
    {
        cout << "ATOM IN MAXUNK<<"<<"\n";
        AccVal=ResAccess[mole->getUnknownRes().getNum()];
        ofstr.setf(ios::right);  ofstr.width(6);       ofstr<<AccVal<<" ";
        ofstr.unsetf(ios::right);ofstr.setf(ios::left);
        oss.str("");
        oss<<"UNK"<<mole->getUnknownRes().getFNum();
        ofstr<< oss.str();
        for (unsigned int N=oss.str().length(); N <= 8; N++) ofstr<<" ";
        ofstr.width(4); ofstr.unsetf(ios::left);ofstr.setf(ios::right);
        for (ItCAtom itA = mole->getUnknownRes().firstAtom();itA != mole->getUnknownRes().lastAtom();++itA)
            {
                if (AtmAccess[(*itA)->getNum()] != 0)
                    {ofstr<< AtmAccess[(*itA)->getNum()]<<" ";break;}
            }
        ofstr<<"GROUP   ";
      ofstr<<"\n";




    }*/
    /*
    typedef struct
    {
        AtomList list;
        std::string name;
        std::string type;
        std::string subtype;
    } AtomSet;
    typedef struct
    {
        BondList list;
        std::string name;
        std::string type;
        std::string subtype;
    } BondSet;
    typedef struct
    {
       ResiduList list;
        std::string name;
        std::string type;
        std::string subtype;
    } ResSet;
*/


    if (!withSets && !withCycles&& NSet!=0) {ofstr<<endl; return;}
    ofstr<< "@<TRIPOS>SET"<<endl;

    if (withSets)
    {
        // Starting with Atom Sets :
        for(std::vector<AtomSet>::iterator itAS =atomsets.begin(); itAS != atomsets.end(); itAS++)
        {
            ofstr<< (*itAS).name<<" STATIC     ATOMS     "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique((*itAS).list);
            for (ItCAtom itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
            {
                const Atom& atm=**itAA;

                if (find(atomlist.begin(),atomlist.end(),&atm) == atomlist.end())
                {

                    continue;
                }
                const size_t pos= std::distance(atomlist.begin(),find(atomlist.begin(),atomlist.end(),&atm))+1;
                //cout <<fname<<"\tATM\t"<< atm.getIdentifier()<<" " << atm.getNum()<< " " << AtmAccess[atm.getNum()]<<endl;
                if (pos==0)continue;
                oss<< pos<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
            }
            ofstr << countTOT<<" "<<oss.str()<<endl;
            SetId++;
        }

    }
    /*for(std::vector<ResSet>::iterator itAS = mole->firstResSet(); itAS != mole->lastResSet(); itAS++)
          {
            ofstr<< (*itAS).name<<" STATIC     SUBSTS    "<<(*itAS).subtype<<endl;
            unsigned int count=1,countTOT=0;
            oss.str("");
            sortAndUnique((*itAS).list);
            for (ItRes itAA = (*itAS).list.begin(); itAA != (*itAS).list.end(); itAA++)
              {
                Residu &res = **itAA;
                if (ResAccess[res.getNum()]==0)continue;
                if (AtmAccess[res.getRootAtom()->getNum()] == 0)continue;
                cout << fname<<"\tRES\t"<<res.getRootAtom()->getIdentifier()<<" " << res.getRootAtom()->getNum()<< " " << AtmAccess[res.getRootAtom()->getNum()]<<endl;
                oss<< AtmAccess[res.getRootAtom()->getNum()]<<" "; count++;countTOT++;
                if (count == 18 && countTOT != (*itAS).list.size()) {oss <<"\\"<<endl;count=0;}
              }
            ofstr << countTOT<<" "<<oss.str()<<endl;
            SetId++;
          }*/


    //if (ossCol.str().length()) ofstr <<ossCol.str();

    // }


}

void MoleWriter::putDataInFile()
{

}
