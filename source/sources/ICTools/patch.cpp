#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICCalcs/fragments.h"
#include "headers/ICMole/bond.h"
#include "headers/ICMole/atom.h"
using namespace std;
using namespace ICMole;

struct Mutants
{
    string RES;
    int  RES_ID;
    string CHAIN;

};

struct Bornes
{
    string ChainName;
    int ATOM_MIN,
    ATOM_MAX,
    UNIP_MIN,
    UNIP_MAX;
};

double calcHScoreO(const Coords& H, const Coords& AtmCD, const Coords &AtmCA, const Coords &AtmD_L, const double& dist, const bool &water)
{
    const double distH_D = H.calcDist(AtmCD,1.2);
    const double distH_A = H.calcDist(AtmCA,dist+0.1);


    const double angl_HDDl = (AtmCD.calcAngle(AtmD_L,H));
    //cout << " ||-> "<<distH_D<< " " << distH_A<<" " <<minAngl<<"<"<< angl_HDDl<<"<"<<maxAngl<<endl;

    //if (angl_HDDl >= maxAngl || angl_HDDl <= minAngl) { return 0;}
    const double anglAHD = (H.calcAngle(AtmCA,AtmCD))*180/PI;              if (anglAHD <=130) { return 0;}
    const double alpha = 34 - distH_A*sin( distH_D*sin(angl_HDDl)/dist);
    const double angl_ADDl=((AtmCD.calcAngle(AtmD_L,AtmCA))*180/PI); if (angl_ADDl > 109.5+alpha || angl_ADDl < 109.5-alpha) {// cout <<angl_ADDl<< " avoid"<<endl;
        return 0;}
    const double diffAHD=180-anglAHD;
    if (diffAHD <= 30) {if (water) return 1+(diffAHD)/50; else return 1;}
    else if (diffAHD <= 80) return 1-(diffAHD-30)/50;
    else return 0;

}


void IChemSwitch::helpPatch() const
{

}
void IChemSwitch::patch() const   throw(ICMole::MoleExcept)
{
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize !=10) {throw MoleExcept(9010301,"IChem::patch","Not enough parameters");}

    const string origProt = Input_Values.at(0);
    const string inProt   = Input_Values.at(1);
    const string inLig    = Input_Values.at(2);
    const string outProt  = Input_Values.at(3);
    const string outLig   = Input_Values.at(4);
    const string outSite  = Input_Values.at(5);
    const string outUniP  = Input_Values.at(6);
    const string outUniS  = Input_Values.at(7);
    const string outInts  = Input_Values.at(8);
    const string NAME     = Input_Values.at(9);
    string Bornestr = "";
    string Mutantstr= "";

    for (std::map<std::string,std::vector<std::string> >::const_iterator
         it = Opt_Values.begin();
         it != Opt_Values.end();
         it++)
    {
        const std::string &opt_name = (*it).first;
        const vector<string> & opt_val = (*it).second;
        const std::string& value = opt_val.at(0);
        //  -rn -cn --values -sim -outInt -match -max -size --all_cliques -score
        if (opt_name.compare("-mutant")            == 0) { Mutantstr=value;}
        else if (opt_name.compare("-bornes")       == 0) { Bornestr=value;}

    }

    vector<Bornes> LISTBORNES;
    vector<Mutants> LISTMUTANT;

    vector<string> tokBornes,tokMutant;
    tokenStr(Bornestr,tokBornes,"-");
    for (size_t i=0; i<tokBornes.size();++i)
    {
        vector<string> tok2;
        tokenStr(tokBornes.at(i),tok2,"_");
        if (tok2.size() != 5) throw MoleExcept(9010401,"IChem::patch","Wrong borne format");
        Bornes BS;
        BS.ChainName=tok2.at(0);
        BS.ATOM_MIN=atoi(tok2.at(1).c_str());
        BS.ATOM_MAX=atoi(tok2.at(2).c_str());
        BS.UNIP_MIN=atoi(tok2.at(3).c_str());
        BS.UNIP_MAX=atoi(tok2.at(4).c_str());
        LISTBORNES.push_back(BS);
    }


    tokenStr(Mutantstr,tokMutant,"-");
    for (size_t i=0; i<tokMutant.size();++i)
    {
        vector<string> tok2;

        tokenStr(tokMutant.at(i),tok2,"_");
        if (tok2.size() != 3) throw MoleExcept(9010401,"IChem::patch","Wrong mutant format");
        Mutants MU;
        MU.CHAIN=tok2.at(0);
        MU.RES = tok2.at(1);
        MU.RES_ID=atoi(tok2.at(2).c_str());
        LISTMUTANT.push_back(MU);
    }

    Complex complex;



    MoleReader mread(inProt);
    mread.loadInComplex(complex,MoleType::PROTEIN);

    mread.loadNewFile(inLig);
    mread.loadInComplex(complex,MoleType::LIGAND);

    Molecule origMole;
    mread.loadNewFile(origProt);
    mread.loadNextMolecule(origMole,MoleType::PROTEIN);

    Molecule &protein=*complex.getMole(MoleType::PROTEIN);
    Molecule &ligand =*complex.getMole(MoleType::LIGAND);


    protein.setName(NAME+"_protein");
    ligand.setName(NAME+"_ligand");

    vector<string> chains;
    map<Residu*,string> restochain;

    string name;
    // Color Ligand atoms :
    // Correct ligand names :
    for (ItCAtom itA = ligand.firstAtom();itA != ligand.lastAtom();++itA)
    {
        Atom& atom =**itA;
        name="";
        for (size_t i=0; i< atom.getName().length();i++)
            if (atom.getName().at(i)!='_')name+=atom.getName().at(i);
        atom.setName(name);
        atom.setColor(212,212,212);
    }

    ResiduList todel,siteList;
    unsigned int NWATER=0;
    for (size_t iRes =0; iRes < protein.numResidus();++iRes)
    {
        Residu& residu = protein.getResidu(iRes);
        for (size_t iAtm=0;iAtm<residu.numAtom();++iAtm)
        {
            Atom& atom =residu.getAtom(iAtm);
            name="";
            for (size_t i=0; i< atom.getName().length();i++)
                if (atom.getName().at(i)!='_')name+=atom.getName().at(i);
            atom.setName(name);
            atom.assignBits();
        }
        // reset chains to missing ones :
        if (residu.getChainName()=="*")
        {
            Atom &firstAtm=residu.getAtom(0);
            for (size_t iAtm=0;iAtm < origMole.numAtom();++iAtm)
            {
                const Atom& comp = origMole.getAtom(iAtm);
                if (firstAtm.fixpos.calcDist(comp.fixpos,0.2)<0.05)
                { cout << firstAtm.getIdentifier()<<"\t"<<comp.getIdentifier()<<endl;
                    chains.push_back(comp.getResidu()->getChainName());
                    restochain.insert(pair<Residu*,string> (&residu,comp.getResidu()->getChainName()));
                }
            }
        }
        if (residu.getName()=="HOH")
        {
            unsigned int NHBond=0;
            Atom* Ox=(Atom*)NULL,*H1=(Atom*)NULL,*H2=(Atom*)NULL;
            for (size_t iAtm=0;iAtm<residu.numAtom();++iAtm)
            {
                Atom& atom =residu.getAtom(iAtm);
                if (atom.isOxygen()) Ox=&atom;
                else if (atom.isHydrogen())
                {
                    if (H1==(Atom*)NULL)H1=&atom; else H2=&atom;
                }
            }

            for (ItCAtom itA = protein.firstAtom();itA != protein.lastAtom();++itA)
            {
                Atom& atm = **itA;
                if (atm.getResidu() == &residu)continue;
                if (!atm.props.isAcceptor() && !atm.props.isDonor())continue;
                if (atm.fixpos.calcDist(Ox->fixpos,3.6)>3.5)continue;
                if (atm.props.isAcceptor())
                {
                    double value= calcHScoreO(H1->fixpos, Ox->fixpos, atm.fixpos, H2->fixpos, atm.fixpos.calcDist(Ox->fixpos,3.6), true);
                    if (value >0) NHBond++;
                    value= calcHScoreO(H2->fixpos, Ox->fixpos, atm.fixpos, H1->fixpos, atm.fixpos.calcDist(Ox->fixpos,3.6), true);
                    if (value >0) NHBond++;
                }
                else if (atm.props.isDonor())
                {
                    Atom *DL=(Atom*)NULL;
                    for (unsigned int i=0; i < atm.getNumBond();i++)
                    {
                        if (!atm.getAtomLinked(i).isHydrogen()){ DL=&atm.getBond(i)->getOtherAtom(atm);break;}
                    }
                    for (unsigned int i=0; i < atm.getNumBond();i++)
                    {
                        if (!atm.getAtomLinked(i).isHydrogen())continue;
                        double value= calcHScoreO(atm.getAtomLinked(i).fixpos, atm.fixpos, Ox->fixpos, DL->fixpos, atm.fixpos.calcDist(Ox->fixpos,3.6), true);
                        if (value >0) NHBond++;
                    }

                }

            }


            if (NHBond < 2){ todel.push_back(&residu);continue;}
            else NWATER++;
        }
        protein.ringPerception();


        double bestdist=1000;
        for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
        {
            const Atom& atmP=residu.getAtom(iPAtm);
            if (atmP.isHydrogen() || atmP.getMOL2Type()=="Du")continue;
            for (size_t iLAtm=0;iLAtm < ligand.numAtom();++iLAtm)
            {
                const Atom& atmL=ligand.getAtom(iLAtm);
                if (atmL.isHydrogen() || atmL.getMOL2Type()=="Du")continue;
                const double dist=atmP.fixpos.calcDist(atmL.fixpos,bestdist+0.1);
                if (dist < bestdist)bestdist=dist;
            }
        }


        // Checking residu distance from ligand :
        if (bestdist >= 6.5)
        {
            for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
            {
                Atom& atmP=residu.getAtom(iPAtm);
                if (atmP.isCarbon())
                {
                    protein.addAtomColor(&atmP);
                    atmP.setColor(255,0,255);
                }
            }
            continue;
        }
        else
        {
            siteList.push_back(&residu);
            // If in site : coloring according to
            // Checking mutant :
            bool isMutant=false;
            if (!LISTMUTANT.empty())
                for (vector<Mutants>::iterator it=LISTMUTANT.begin();it != LISTMUTANT.end();it++)
                {
                    if (residu.getChainName() != (*it).CHAIN)continue;
                    if (residu.getName() != (*it).RES) continue;
                    if (residu.getFNum()== (*it).RES_ID) isMutant=true;
                }
            if (isMutant) cout << endl<<NAME<<"|IS MUTANT"<<endl;
            for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
            {
                Atom& atmP=residu.getAtom(iPAtm);
                if (atmP.isCarbon())
                {
                    protein.addAtomColor(&atmP);
                    if (!isMutant) atmP.setColor(0,255,0);
                    else atmP.setColor(255,255,  0);
                }
            }

        }
    }


    for (ItRes it = todel.begin(); it != todel.end(); it++)
    {
        protein.delResidu(**it);

    }
    cout<< NAME<<"|NWAT\t"<<NWATER<<endl;

    sort(chains.begin(),chains.end());
    vector<string>::iterator it=unique(chains.begin(),chains.end());
    chains.resize(std::distance(chains.begin(),it));

    if (chains.size()==1)
    {
        Residu &res = *(*restochain.begin()).first;
        res.getChain()->setName(chains.at(0));
    }

    for (ItCBond itCB = protein.firstBond();itCB != protein.lastBond();++itCB)
    {
        Bond &bd= **itCB;
        bd.assignBits();
    }
    MoleWriter mw(outProt);mw.writeMOL2(&protein);
    mw.newFile(outLig);mw.writeMOL2(&ligand);



    for (size_t iRes=0;iRes< protein.numResidus();++iRes)
    {
        Residu& residu=protein.getResidu(iRes);
        if (find(siteList.begin(),siteList.end(),&residu)==siteList.end()) residu.setUse(false);
        else
        {
            residu.setUse(true);
            for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
            {
                Atom& atmP=residu.getAtom(iPAtm);
                if (!atmP.isMetallic())continue;
                if (atmP.getResidu()->numAtom()==1 &&
                        atmP.getName().length() < atmP.getResiduName().length())
                {
                    atmP.getResidu()->setFNum(atoi(atmP.getResiduName().substr(atmP.getName().length()).c_str()));
                    atmP.getResidu()->setName(atmP.getResiduName().substr(0,atmP.getName().length()));
                }
            }
        }
    }

    protein.setName(NAME+"_site");
    mw.newFile(outSite);
    mw.setOnlyUsed(true);
    mw.writeMOL2(&protein);
    Interactions ints(complex);
    InterResults intsresults;
    ints.calcInteractions(ligand,intsresults,true);
    ints.interToMOL2(intsresults,false,false,false,true);
    mw.newFile(outInts+"_M.mol2");
    intsresults.Ints.setName(NAME+"_INTS_M");
    unsigned short CA,CZ,O,OD1,N,NZ,Zn;
    CA=0;CZ=0;O=0;OD1=0;N=0;NZ=0;Zn=0;
    cout << intsresults.Ints.numAtom()<<endl;
    for (size_t iAtm=0;iAtm < intsresults.Ints.numAtom();++iAtm)
    {
        const Atom& atm =  intsresults.Ints.getAtom(iAtm);

        if (atm.getName()=="CA")CA++;
        else if (atm.getName()=="CZ")CZ++;
        else if (atm.getName()=="O")O++;
        else if (atm.getName()=="OD1")OD1++;
        else if (atm.getName()=="N")N++;
        else if (atm.getName()=="NZ")NZ++;
        else if (atm.getName()=="Zn")Zn++;

    }
    cout << NAME<<"|INTS\tMERGE\t"
         << (CA+CZ+O+OD1+N+NZ+Zn)<<"\t"
         << CA<<"\t"
         << CZ<<"\t"
         << O<<"\t"
         << OD1<<"\t"
         << N<<"\t"
         << NZ<<"\t"
         << Zn<<"\n";
    mw.writeMOL2(&intsresults.Ints);
    ints.interToMOL2(intsresults,false,false,true,false);
    intsresults.Ints.setName(NAME+"_INTS_C");
    CA=0;CZ=0;O=0;OD1=0;N=0;NZ=0;Zn=0;
    for (size_t iAtm=0;iAtm < intsresults.Ints.numAtom();++iAtm)
    {
        const Atom& atm =  intsresults.Ints.getAtom(iAtm);
        if (atm.getName()=="CA")CA++;
        else if (atm.getName()=="CZ")CZ++;
        else if (atm.getName()=="O")O++;
        else if (atm.getName()=="OD1")OD1++;
        else if (atm.getName()=="N")N++;
        else if (atm.getName()=="NZ")NZ++;
        else if (atm.getName()=="Zn")Zn++;

    }
    cout << NAME<<"|INTS\tCENTER\t"
         << (CA+CZ+O+OD1+N+NZ+Zn)<<"\t"
         << CA<<"\t"
         << CZ<<"\t"
         << O<<"\t"
         << OD1<<"\t"
         << N<<"\t"
         << NZ<<"\t"
         << Zn<<"\n";
    mw.newFile(outInts+"_C.mol2");
    mw.writeMOL2(&intsresults.Ints);
    ints.interToMOL2(intsresults,true,false,false,false);
    intsresults.Ints.setName(NAME+"_INTS_L");
    mw.newFile(outInts+"_L.mol2");
    mw.writeMOL2(&intsresults.Ints);
    CA=0;CZ=0;O=0;OD1=0;N=0;NZ=0;Zn=0;
    for (size_t iAtm=0;iAtm < intsresults.Ints.numAtom();++iAtm)
    {
        const Atom& atm =  intsresults.Ints.getAtom(iAtm);
        if (atm.getName()=="CA")CA++;
        else if (atm.getName()=="CZ")CZ++;
        else if (atm.getName()=="O")O++;
        else if (atm.getName()=="OD1")OD1++;
        else if (atm.getName()=="N")N++;
        else if (atm.getName()=="NZ")NZ++;
        else if (atm.getName()=="Zn")Zn++;

    }
    cout << NAME<<"|INTS\tLIG\t"
         << (CA+CZ+O+OD1+N+NZ+Zn)<<"\t"
         << CA<<"\t"
         << CZ<<"\t"
         << O<<"\t"
         << OD1<<"\t"
         << N<<"\t"
         << NZ<<"\t"
         << Zn<<"\n";
    ints.interToMOL2(intsresults,false,true,false,false);
    intsresults.Ints.setName(NAME+"_INTS_P");
    CA=0;CZ=0;O=0;OD1=0;N=0;NZ=0;Zn=0;
    for (size_t iAtm=0;iAtm < intsresults.Ints.numAtom();++iAtm)
    {
        const Atom& atm =  intsresults.Ints.getAtom(iAtm);
        if (atm.getName()=="CA")CA++;
        else if (atm.getName()=="CZ")CZ++;
        else if (atm.getName()=="O")O++;
        else if (atm.getName()=="OD1")OD1++;
        else if (atm.getName()=="N")N++;
        else if (atm.getName()=="NZ")NZ++;
        else if (atm.getName()=="Zn")Zn++;

    }
    cout << NAME<<"|INTS\tPROT\t"
         << (CA+CZ+O+OD1+N+NZ+Zn)<<"\t"
         << CA<<"\t"
         << CZ<<"\t"
         << O<<"\t"
         << OD1<<"\t"
         << N<<"\t"
         << NZ<<"\t"
         << Zn<<"\n";
    mw.newFile(outInts+"_P.mol2");
    mw.writeMOL2(&intsresults.Ints);
    ofstream ofs;ofs.open((outInts+"_LOG").c_str());
    ofs<<ints.toString(intsresults);
    ofs.close();
    ints.genIFP(intsresults);
    ofs.open((NAME+"_IFP").c_str());
    ofs<<intsresults.IFPString<<"\n";
    ofs<<intsresults.IFP.toString()<<"\n";
    ofs.close();
    cout << NAME<<"|IFP\tHEADER\t"<<intsresults.IFPString<<endl;
    cout << NAME<<"|IFP\tVALUE\t"<<intsresults.IFP.toString()<<endl;

    FragGenerator FragGen(complex);
    FragGen.fragmentLigand();
    FragGen.calcInteractions();
    FragGen.saveFragments(NAME);

    if (LISTBORNES.empty()) return ;



    protein.setName(NAME+"_protein");

    todel.clear();

    for (size_t iRes=0;iRes< protein.numResidus();++iRes)
    {
        Residu& residu = protein.getResidu(iRes);
        const    int& resNum    = residu.getFNum();
        const string& chainName = residu.getChainName();
        if (residu.getResType()!=ResType::STD_AA && residu.getResType() != ResType::MOD_AA)continue;
        // Finding the rule :
        vector<Bornes>::const_iterator GoodBorne= LISTBORNES.end();
        for (vector<Bornes>::const_iterator itBorne = LISTBORNES.begin(); itBorne != LISTBORNES.end(); itBorne++)
        {
            const Bornes& testBorne = *itBorne;
            if (testBorne.ChainName != chainName)continue;
            // cout << chainName << resNum<< " " << testBorne.ATOM_MIN<< " " << testBorne.ATOM_MAX<<endl;
            if (resNum < testBorne.ATOM_MIN) continue;
            if (resNum > testBorne.ATOM_MAX+1)continue;
            GoodBorne=itBorne;
        }
        if (residu.getResType() == ResType::STD_AA && GoodBorne == LISTBORNES.end())
        { todel.push_back(&residu);
            cout << NAME<<"|UNIPROT|PROTEIN|DELETION\t" <<residu.getIdentifier()<<endl;
            continue;
        }
        residu.setFNum((*GoodBorne).UNIP_MIN+(residu.getNum()-(*GoodBorne).ATOM_MIN));

    }
    for (ItRes it = todel.begin(); it!= todel.end();it++)
    {
        protein.delResidu(**it);
    }

    protein.clearSets();
    for (size_t iRes =0; iRes < protein.numResidus();++iRes)
    {
        Residu& residu = protein.getResidu(iRes);


        double bestdist=1000;
        for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
        {
            const Atom& atmP=residu.getAtom(iPAtm);
            if (atmP.isHydrogen() || atmP.getMOL2Type()=="Du")continue;
            for (size_t iLAtm=0;iLAtm < ligand.numAtom();++iLAtm)
            {
                const Atom& atmL=ligand.getAtom(iLAtm);
                if (atmL.isHydrogen() || atmL.getMOL2Type()=="Du")continue;
                const double dist=atmP.fixpos.calcDist(atmL.fixpos,bestdist+0.1);
                if (dist < bestdist)bestdist=dist;
            }
        }


        // Checking residu distance from ligand :
        if (bestdist >= 6.5)
        {
            for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
            {
                Atom& atmP=residu.getAtom(iPAtm);
                if (atmP.isCarbon())
                {
                    protein.addAtomColor(&atmP);
                    atmP.setColor(255,0,255);
                }
            }
            continue;
        }
        else
        {
            // If in site : coloring according to
            // Checking mutant :
            bool isMutant=false;
            if (!LISTMUTANT.empty())
                for (vector<Mutants>::iterator it=LISTMUTANT.begin();it != LISTMUTANT.end();it++)
                {
                    if (residu.getChainName() != (*it).CHAIN)continue;
                    if (residu.getName() != (*it).RES) continue;
                    if (residu.getFNum()== (*it).RES_ID) isMutant=true;
                }
            if (isMutant) cout << endl<<NAME<<"|IS MUTANT"<<endl;
            for (size_t iPAtm=0;iPAtm < residu.numAtom();++iPAtm)
            {
                Atom& atmP=residu.getAtom(iPAtm);
                if (atmP.isCarbon())
                {
                    protein.addAtomColor(&atmP);
                    if (!isMutant) atmP.setColor(0,255,0);
                    else atmP.setColor(255,255,  0);
                }
            }

        }
    }
    mw.setOnlyUsed(false);
    mw.newFile(outUniP);mw.writeMOL2(&protein);
    // Generating Protein based on Uniprot :
    mw.setOnlyUsed(true);
    mw.newFile(outUniS);mw.writeMOL2(&protein);

}

