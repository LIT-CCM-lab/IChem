//#define ICHEM_DEBUG
#include "headers/ICCalcs/ppi.h"
#include "headers/ICTools/switch.h"
#include "headers/ICMole/grid.h"
#include "headers/ICMole/box.h"
#include "headers/ICCalcs/volsite.h"

using namespace std;
using namespace ICMole;



unsigned int verbose = 1;
bool w_chain=false,w2_chain=false;
std::vector<std::string> vecChains;
int nchain=0;
string ChainSele,ChainSele2;
int nChaintot =1;
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/*!<  \fn PPI::PPI(vector<Chain> vecChain)
 * \brief PPI constructor
 *
 * \param vecVhain : Chain of the interface
*/

PPI::PPI(Complex &cp, string name):cp_orig(cp),PDB_name(name)
{
#ifdef ICHEM_DEBUG
    cout << "    |DETPPI|START\t"<<endl<<PDB_name<<"|DETPPI|HEADER\t##### detectPPi  "<<endl;
#endif
    /*vecChain = chains;
    resInter = resInt;*/

}





//    Molecule mole,mole2 ;
//    ResiduList reslist;
//    std::copy(mole.getChain(0).first(),mole.getChain(0).last(),reslist.end());
//    mole.moveResidu(reslist,mole2);



//bool comparaisoninterface (InterResults* i,InterResults* j) {
//    return (i->Ints.getMaxNumAtom() >j->Ints.getMaxNumAtom());
//}

bool comparaisoninterface ( PPI::Interface *i, PPI::Interface *j){
    return (i->size > j->size);
}


void PPI::mapinitialisation(){

    // BOUCLE 1~ERE CHAINE
    //% verbose a faire initialisation
    int chaintreat=0;
#ifdef ICHEM_DEBUG
    cout << PDB_name<<"|DETPPI|HEADER\t##### Map Initialisation  "<<endl;
    cout << PDB_name<<"|DETPPI|HEADER\t Nombre de chaine données ; "<< nchain << endl;
    //% A revoir si on donne l'unite biologique
    if (cp_orig.getNumMolecule() > 15 && !w_chain && !w2_chain){
        cout << PDB_name<<"|DETPPI|HEADER\t Too much Chains in the PDB" << endl;
        //return;
    }
#endif
    for (ItCMole itMole = cp_orig.firstMole() ; itMole!=cp_orig.lastMole();++itMole)
    {
        const Molecule& mole = **itMole;
        if ((nchain==1 && mole.getName()!= vecChains[0]) || (nchain>1 && find(vecChains.begin(),vecChains.end(),mole.getName())==vecChains.end())) {
            continue;
        }
        map<string,vector<PPI::Interaction> > mapvector;
        for(ItCMole itMole2=cp_orig.firstMole() ; itMole2!=cp_orig.lastMole();++itMole2){
            const Molecule& mole2 = **itMole2;
            if (mole2.getName()==mole.getName()) continue;

            //% TRouver un moyen de vérifier la longueru des chaines ...
            if ((nchain>1 && find(vecChains.begin(),vecChains.end(),mole2.getName())==vecChains.end()) ) continue;
#ifdef ICHEM_DEBUG
            cout << PDB_name<<"|DETPPI|MAPINI\tBetween \t" << mole.getName() << " || " << mole2.getName() << endl ;

#endif
            vector<PPI::Interaction> T;
            mapvector.insert(pair<string,vector<PPI::Interaction> >(mole2.getName(),T));
        }
        // ajout des vecteurs receptacle d'interaction dans la map
        mapinteraction.insert(pair<string,map<string,vector<PPI::Interaction> > >(mole.getName(),mapvector));
        chaintreat++;
//        mole.setUse(false);
    }

}

void PPI::mapdetection(){


    // Vu que l'on ne parcours plus la seconde chaine il se peux que celle ci soit désactiver plus tot
    Grid &bigGrid=cp_orig.genGrid(4.5);
    AtomList adjacents;
    vector<string> chaintreat;
    for (ItCMole itMole = cp_orig.firstMole() ; itMole!=cp_orig.lastMole();++itMole)
    {
        Molecule& mole = **itMole;

        //        cout << "mole : "<< mole.getName() << endl;

        if ((nchain==1 && mole.getName()!= vecChains[0]) || (nchain>1 && find(vecChains.begin(),vecChains.end(),mole.getName())==vecChains.end())) {
            mole.setUse(false); //% Voir quoi faire ici (decommenter le 26/02/2015
            continue;
        }
        map<string,vector<PPI::Interaction> > mapvector;

        string chainName = mole.getName();
        for (ItCAtom itAtom = mole.firstAtom();  itAtom!=mole.lastAtom(); ++itAtom ){
            Atom& atom = **itAtom;
            bigGrid.getAdjacentAtoms(adjacents,atom,5,true);
            //            cout << "liste adjacences: " << adjacents.size() << endl;
            for (ItCAtom itAtombox = adjacents.begin();itAtombox != adjacents.end() ;++itAtombox ){
                Atom& atomCheck= **itAtombox;
                vector<PPI::Interaction> T;
                if (atomCheck.getParent().getName() != chainName ){
#ifdef ICHEM_DEBUG

                    cout << chainName << " vs " << atomCheck.getParent().getName() << endl;
                    cout << "   First Atom : "<< atom.getIdentifier()<< "   ||    " <<atomCheck.getIdentifier() << " dist :" << atomCheck.fixpos.calcDist(atom.fixpos) << endl;
#endif
                    Interaction inter;
                    inter.atomA= &atom;
                    inter.atomB= &atomCheck;
                    //                    if (inter.atomA->getName()=="DuAr" && find(vecChains.begin(),vecChains.end(),(inter.atomA->getParent().getName())) == vecChains.end()) {inter.atomA->getResidu()->setUse(false);}
                    //                    if (inter.atomB->getName()=="DuAr" && find(vecChains.begin(),vecChains.end(),(inter.atomB->getParent().getName())) == vecChains.end()) {inter.atomB->getResidu()->setUse(false);}
                    inter.dist=atomCheck.fixpos.calcDist(atom.fixpos);
                    //                    if (find(chaintreat.begin(),chaintreat.end(),atomCheck.getResidu()->getChainName())!=chaintreat.end()) cout << "j'accepte pas" << endl; continue;
                    mapinteraction[atom.getParent().getName()][atomCheck.getParent().getName()].push_back(inter);
                    //                    cout << inter.atomA->getParent().getName() << " avec " <<inter.atomB->getParent().getName()<<endl;
                    //                    cout << " ajout a : " << atom.getParent().getName()<<" - " << atomCheck.getParent().getName() << endl;
                    //                    cout << "nombre d'interaction " << mapinteraction[atom.getParent().getName()][atomCheck.getParent().getName()].size()<< endl;
                    //                    if (atomCheck.getName()=="H") cout << "ajout hydrogen" << endl;
                    atomCheck.getParent().setUse(false);
                }
            }
            atom.getParent().setUse(false);
        }
        //        chaintreat.push_back(mole.getName());
        //        cout << "size " << mapinteraction["H"]["N"].size() << endl;
        mole.setUse(false);
    }

}

void PPI::makeinterface(){
#ifdef ICHEM_DEBUG
    cout << PDB_name<<"|DETPPI|HEADER\t##### Make Interface  "<<endl;
#endif
    /// boucle de creation des objets interfaces
    ///  fonction qui ressort vector d'interface non?
    ///


#ifdef ICHEM_DEBUG
    cout << PDB_name<<"|DETPPI|HEADER\t##### Creation des interfaces detectées (+Intergene)"<<endl;
#endif


    vector<string> chaintreat;
    for (ItCMole itMole = cp_orig.firstMole() ; itMole!=cp_orig.lastMole();++itMole)
    {

        Molecule& chain1 = **itMole;

//        cout << chain1.toString() << endl;
        if ((nchain==1 && chain1.getName()!= vecChains[0]) || (nchain>1 && find(vecChains.begin(),vecChains.end(),chain1.getName())==vecChains.end())) {continue;}
        const int taille=10;
        if (chain1.getName()=="XX" ){continue;}
        //        if (chain1.getCountRes() < taille ){ continue;}
        //        ItCMole itMole2(itMole);
        for (ItCMole itMole2=cp_orig.firstMole();itMole2!=cp_orig.lastMole();++itMole2)
        {
            Molecule& chain2 = **itMole2;

//            cout << endl<< endl << "NEw Chain" << endl << chain2.toString() << endl;

            if (find(chaintreat.begin(),chaintreat.end(),chain2.getName())!=chaintreat.end()) continue;
            if (chain2.getName()=="XX" or chain1.getName()==chain2.getName())continue;
            if ((nchain>1 && find(vecChains.begin(),vecChains.end(),chain2.getName())==vecChains.end()) ) { //% trouver un moyen de verifier la taille
#ifdef ICHEM_DEBUG
                cout <<"on oublie car a "<< chain2.getName() << " en tant que chain2 " << endl;
#endif
                continue;
            }

            // Verbose
            if ((mapinteraction.at(chain1.getName()).at(chain2.getName()).size()<20) && (mapinteraction.at(chain2.getName()).at(chain1.getName()).size()<20)){
                cout << "Interaction trop petite entre " << chain1.getName() << " et " << chain2.getName() << endl;
                continue;
            }

            AtomList toComp;

            for( std::vector<Interaction>::iterator i = mapinteraction.at(chain1.getName()).at(chain2.getName()).begin();i!=mapinteraction.at(chain1.getName()).at(chain2.getName()).end();++i)
            {
                toComp.push_back( (*i).atomA);
                toComp.push_back( (*i).atomB);
                (*i).atomB->getResidu()->setUse(true);
                (*i).atomA->getResidu()->setUse(true); // commenté dans PIC
//                cout << (*i).atomA->getResidu()->getIdentifier() << endl;
//                cout << (*i).atomB->getResidu()->getIdentifier() << endl;
            }
            if (toComp.size()==0){
                cout <<"Pas d'atome en interaction" << endl;
                continue;
            }

            //// tri et filtre liste d'atome (peut etre pas necessaire
            SortUniqueVector(toComp);

            //// recupere les residus des atomes de la liste faire une fonction transforme vecteur atom en vecteur residu
#ifdef ICHEM_DEBUG
            cout << PDB_name<<"|DETPPI|MAKINT|D\t Passage des atomes aux residus "<< chain1.getName() << " || "<< chain2.getName() <<endl;
            for (ItAtom atomV = toComp.begin() ;atomV != toComp.end();atomV++){
                cout << PDB_name<<"|DETPPI|MAKINT|D\t"
                     << "Res Call :\t" << (*atomV)->getResidu()->getName() << "  "<< (*atomV)->getResidu()->getNum() << "  " << (*atomV)->getResidu()->getChainName()<< " // " <<  endl;
            }
#endif


            MoleWriter mwtest;
            string A=chain1.getName();
            string B=chain2.getName();
//            mwtest.newFile(PDB_name+"_site_"+A+"_"+B+".mol2");
//            mwtest.setOnlyUsed(true);

//            mwtest.writeMOL2(cp_orig,PDB_name);

            chain1.setUse(true);
            chain2.setUse(true);

            chain1.ringPerception(); // On ne peut donc pas calculer toutes les interfaces d'un coup
            chain2.ringPerception();
            MoleWriter mw;
            Interactions ins(cp_orig);
            InterResults results;
//            string A=chain1.getName();
//            string B=chain2.getName();
            mw.newFile(PDB_name+"_prot_"+A+"_"+B+".mol2");
            mw.setOnlyUsed(true);


            mw.writeMOL2(cp_orig,PDB_name);


            ins.getGrid().setWHydrogen(true);
//            chain2.setMoleType(MoleType::LIGAND);
            ins.setDist_Arom(4.5);
            ins.calcInteractionsppi(results,true,true);
            //        cout << ins.toString(results);
//            chain2.setMoleType(MoleType::PROTEIN);

            results.chainA=A;
            results.chainB=B;

            if (results.listInters.size() < 5 ){ //% Sortie a formater
                //                    interfaces.push_back(interCheck);
                cout << results.listInters.size() << " interactions for " << PDB_name << " chain : "<< A <<" "<<  B<< endl;
                continue;
            }
            cout << results.listInters.size() << " interactions for " << PDB_name << " chain : "<< A <<" "<<  B<< endl;
            //        cout << ins.toString(results) << endl;



            //            MoleWriter mw;
            ins.interToMOL2(results,false,false,true,false);
            mw.newFile(PDB_name+"_ints_"+A+"_"+B+".mol2");
            results.Ints.setName(PDB_name);
            mw.writeMOL2(&results.Ints);
            //            cout << "valeur de rotpos et fixpos : "<<  results.Ints.getRotpos().toString() << " " << results.Ints.getFixpos().toString() << endl;

            ////Ecris le fichier de log
            string logfile = PDB_name+"_interaction_"+chain1.getName()+"_"+chain2.getName()+".ints"; //*$
            //            string logfile = "interaction.ints"; //nom du fichier d'interaction pour le serveur web

            // a voir si je met un argument pour la sortie des fichiers.
            ofstream of;
            if (fopen(logfile.c_str(),"r")==NULL)
            {
                of.open(logfile.c_str(),ios::out|ios::app);
            }
            else
            {
                of.open(logfile.c_str(),ios::out|ios::app);
            }
            if (of.is_open())
            {
                of<<ins.toString(results);
            }
            of.close();


            Interface* terface= new Interface(chain1.getName(),chain2.getName(),results.listInters.size());

            interfaces.push_back(&ins);
            intersize.push_back(terface);
            scores.push_back(&results);

            // ajout du type d'interaction backbone ou chaine latérale
            int ii= scores.size()-1;
            int nb_inter;
            interfacetype[(ii*3)+0]=0;
            interfacetype[(ii*3)+1]=0;
            interfacetype[(ii*3)+2]=0;

            for (std::vector<InterPoint>::const_iterator it=results.listInters.begin();it != results.listInters.end();it++)
            {
                const InterPoint& interpt=*it;

                if (interpt.interaction == InterType::METAL || interpt.interaction == InterType::PICATION || interpt.interaction == InterType::WHBOND_LIG || interpt.interaction == InterType::WHBOND_PROT || interpt.interaction == InterType::ARFACEFACE || interpt.interaction == InterType::AREDGEFACE) continue;
                else if ( (interpt.Prot_Ref->getName()=="CA" ||interpt.Prot_Ref->getName()=="N" ||interpt.Prot_Ref->getName()=="O") &&  (interpt.Lig_Ref->getName() =="CA" ||interpt.Lig_Ref->getName() =="N" ||interpt.Lig_Ref->getName() =="O")) interfacetype[(ii*3)+0]++;
                else if (!(interpt.Prot_Ref->getName()=="CA" ||interpt.Prot_Ref->getName()=="N" ||interpt.Prot_Ref->getName()=="O") && !(interpt.Lig_Ref->getName() =="CA" ||interpt.Lig_Ref->getName() =="N" ||interpt.Lig_Ref->getName() =="O")) interfacetype[(ii*3)+1]++;
                else interfacetype[(ii*3)+2]++;
                nb_inter++;
            }
#ifdef ICHEM_DEBUG
            cout << "Backbone interactions :" << " B-B : " <<  interfacetype[(ii*3)+0] << " S-S : " << interfacetype[(ii*3)+1] << " B-S : " << interfacetype[(ii*3)+2] << endl;
#endif





            //important
            chain2.setUse(false);


        }
        //        cout <<"score.size() "<< scores.size() << " with chain 1 : "<< chain1.getName() << endl;
        chain1.setUse(false);
        chaintreat.push_back(chain1.getName());

    }

    //% tous déplacer ici





}




void PPI::calcenfouiss(bool outputall) // calcul l'enfouissement des interfaces
{



    for (ItCMole itMole = cp_orig.firstMole() ; itMole!=cp_orig.lastMole();++itMole)
    {
        Molecule& chain1 = **itMole;
        chain1.setUse(true);
    }
    VolSite volsite(cp_orig);
    volsite.proj(55-20,false);
    MoleReader mread2;


    //% Reintergrer le paramètre du choix du nombre d'interface a traité, impecable bien
    for( int i =0;i < intersize.size() ;i++)
    {
        string A = intersize[i]->chain1;
        string B = intersize[i]->chain2;
//        ins.toString(results)
//        cout << "toto "<< interfaces[i]->toString(*scores[i]) << endl;
        ofstream offenf;
        if (outputall){
            string enffile = PDB_name+"_enfouiss_"+A+"_"+B+".enfs"; //*$
            //            string logfile = "interaction.ints"; //nom du fichier d'interaction pour le serveur web

            // a voir si je met un argument pour la sortie des fichiers.

            if (fopen(enffile.c_str(),"r")==NULL)
            {
                offenf.open(enffile.c_str(),ios::out|ios::app);
            }
            else
            {
                offenf.open(enffile.c_str(),ios::out|ios::app);
            }


            offenf<<"enfouissement des interactions" <<endl;
        }

        //        cout << A << "-" <<  B << " avec "  << intersize[i]->size << endl;
        //        cout << PDB_name+"_ints_"+A+"_"+B+".mol2" << endl;
        Molecule test;
        mread2.loadNewFile(PDB_name+"_ints_"+A+"_"+B+".mol2");
        mread2.loadNextMolecule(test,MoleType::LIGAND);


        //        cout << "ma proteine " << test.toString() << endl;
        //        return;

        vector<double> FGPCount;
        vector<vector<double> > FGPBuried;
        for (size_t j=0;j<8;++j)
        {
            FGPCount.push_back(0);
            FGPCount[0] = 0;
            if (j==0)continue;
            vector<double> tmp;
            for (size_t j=0;j<10;++j)tmp.push_back(0); FGPBuried.push_back(tmp);
        }
        for (size_t iAtm=0; iAtm < test.numAtom();++iAtm)
        {
            const Atom& atm= test.getAtom(iAtm);
            Box* box = volsite.getGrid().getBox(atm.fixpos);
            //cout << atm.getName()<<"\t"<<box->getId()<< " " << atm.fixpos.calcDist(box->fixpos)<< " " << volsite.getProjValue(*box)<<endl;
            double pos= floor((volsite.getProjValue(*box)-30)/10);
            if (outputall)
            {
                offenf  << atm.getName()<< " " << volsite.getProjValue(*box)  << endl;
            }

            //            cout << volsite.getProjValue(*box) << "\t col : " << pos+1 << endl;
            //pas de distinction protein-lignadd
            //passer tous les <5 a <7
            if (atm.getName()=="CA")      {FGPBuried[0][pos]+=1;FGPCount[1]++;FGPCount[0]++;}
            else if (atm.getName()=="O")  {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
            else if (atm.getName()=="OD1"){FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
            else if (atm.getName()=="N")  {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
            else if (atm.getName()=="NZ") {FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
            else if (atm.getName()=="CZ") {FGPBuried[3][pos]+=1;FGPCount[4]++;FGPCount[0]++;}
            else if (atm.getName()=="Zn") {FGPBuried[4][pos]+=1;FGPCount[5]++;FGPCount[0]++;}
        }
        for (size_t i=0; i<4;++i)
        {

            for (size_t j=0;j<10;++j)
            {
                if (FGPBuried[i][j]==0)continue;
                //                    cout << i << ":" << j << "valeur " << FGPBuried[i][j]<< endl;
                FGPBuried[i][j]=FGPBuried[i][j]/FGPCount[i+1]*100.0;

            }
            FGPCount[i+1]=FGPCount[i+1]/FGPCount[0]*100.0;

        }
        cout << A << "-" <<B<< "\t";
        //    cout <<cp_orig.ge<<"\t";
        ofstream ofscore;
        ofscore.open("./descriptor.sre", ios::out|ios::app);
        for (size_t i=0;i<=4;++i) {
            cout << FGPCount[i]<< " ";
            ofscore << FGPCount[i]<< " ";
        }
        for (size_t i=0; i<4;++i) // 4 sans metal 5 avec metal
        {

            for (size_t j=0;j<10;++j)
            {

                cout << FGPBuried[i][j]<< " ";
                ofscore <<FGPBuried[i][j]<< " ";
            }

        }
        // ajout du type d'interaction backbone ou chaine latérale
        //        cout << interfacetype[(i*3)+0] << " "<< interfacetype[(i*3)+1] << " "<< interfacetype[(i*3)+2] << " ";
        cout <<"#"<< PDB_name<< A << B <<endl;
        //        ofscore << interfacetype[(i*3)+0] << " "<< interfacetype[(i*3)+1] << " "<< interfacetype[(i*3)+2] << " ";
        ofscore <<"#"<< PDB_name<< "_"<<A  << B << endl;
        ofscore.close();
        if (outputall) offenf.close();

    }



    return;
}

void PPI::detectcavity(bool outputall) // Process
{

    VolSite volsitecav(cp_orig);
    cout << "All chain" << endl;
    volsitecav.setName(PDB_name);
    volsitecav.setMinSize(35);
    volsitecav.proj(55-20,false);

    volsitecav.treatments(5,2.5,false);
//    volsitecav.treatments(5,2.5,outputall);
    volsitecav.clusterCAV(55);
    volsitecav.output_cavs_ppi("all",0);

    for (ItCMole itMole = cp_orig.firstMole() ; itMole!=cp_orig.lastMole();++itMole)
    {
        Molecule& chain1 = **itMole;
        chain1.setUse(false);
    }
//    for (ItCMole itMole = cp_orig.firstMole() ; itMole!=cp_orig.lastMole();++itMole)
//    {

//        Molecule& chain1 = **itMole;
//        chain1.setUse(true);
//        if ((nchain==1 && chain1.getName()!= vecChains[0]) || (nchain>1 && find(vecChains.begin(),vecChains.end(),chain1.getName())==vecChains.end())) {continue;}
//        const int taille=10;
//        if (chain1.getName()=="XX" ){continue;}

//        VolSite volsitecav(cp_orig);
//        cout << chain1.getName() << endl;
//        volsitecav.setName(PDB_name);
//        volsitecav.setMinSize(35);
//        volsitecav.proj(55-20,false);

//        volsitecav.treatments(5,2.5,outputall);
//        volsitecav.clusterCAV(55);
//        string name =chain1.getName();
//        volsitecav.output_cavs_ppi(name,0);
//        chain1.setUse(false);
//    }


//    Grid grid;
//    VolSite volsite(*cp_orig.getMole(MoleType::PROTEIN),grid,1.5);
//    volsite.setName(cp_orig.getMole(MoleType::PROTEIN)->getName());

//    volsite.proj(55-20,false);

//    volsite.treatments(5);
//    volsite.clusterCAV(55);
//    volsite.output_cavs_all("test",0);
    return;
}

bool appartientppi( std::string name, const std::string& liste)
{
    if (name.find(" ,") == string::npos && liste.find(name) != string::npos) return true;
    std::string::size_type lastPos = liste.find_first_not_of(" ,", 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = liste.find_first_of(" ,", lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        if(name == liste.substr(lastPos, pos - lastPos))
            return true;
        // Skip delimiters.  Note the "not_of"
        lastPos = liste.find_first_not_of(" ,", pos);
        // Find next "non-delimiter"
        pos = liste.find_first_of(" ,", lastPos);
    }

    return false;
}




void PPI::defPropsppi(Molecule& protein)
{
    string  resname;
//    cout <<" je suis dans defpropsppi " << protein.getName() << endl;
    for (ItCAtom itA = protein.firstAtom(); itA != protein.lastAtom(); ++itA)
    {
        Atom &atm = **itA;
        const std::string& name = atm.getName();
        PhysProp p = atm.props;
        //
        if (atm.getResidu() != (Residu*)NULL) resname = atm.getResidu()->getName();else resname="";
//        cout << "Atom passé dans defprops : "<< atm.getIdentifier() << "  " << atm.props.toString() << endl ;
        atm.props.clear();
        if (atm.isHydrogen()) continue;
        if (appartientppi(name,"CB") && resname != "CA")                       {
            if (appartientppi(resname,"ALA, VAL, ILE, LEU, PRO, CYS, MET, PHE, TYR, TRP")) {atm.props.setHydrophobic(true);continue;}}
        if (name.compare("CD") == 0 && resname.compare("PRO"))     {atm.props.setHydrophobic(true);continue;}
        if (name.compare("CD1")==0){
            if (appartientppi(resname,"ILE, LEU"))                     {atm.props.setHydrophobic(true);continue;}
            else if (appartientppi(resname, "PHE, TRP, TYR"))                  {atm.props.setAromatic(true);continue;}}
        if (name.compare("CD2")==0){
            if (resname.compare("LEU")==0)                          {atm.props.setHydrophobic(true);continue;}
            else if (appartientppi(resname, "PHE, TRP, TYR"))                  {atm.props.setAromatic(true);continue;}
            else if (appartientppi(resname, "HIS,HID,HIE"))                    {atm.props.setAromatic(true);
                atm.props.setHydrophobic(false);continue;}}
        if (name.compare("CE")==0 && resname.compare("MET"))            {atm.props.setHydrophobic(true);continue;}
        if (name.compare("CE1")==0){
            if (appartientppi(resname,"HIS,HID,HIE"))                      {atm.props.setAromatic(true);
                atm.props.setHydrophobic(false);continue;}
            else if (appartientppi(resname,"PHE, TYR"))                    {atm.props.setAromatic(true);continue;}}
        if (name.compare("CE2")==0 && appartientppi(resname, "PHE, TRP, TYR"))     {atm.props.setAromatic(true);continue;}
        if (name.compare("CE3")==0 && resname.compare("TRP")==0)                {atm.props.setAromatic(true);continue;}
        if (name.compare("CG")==0)
        {    if (appartientppi(resname,"LEU, MET, PRO")){ atm.props.setHydrophobic(true);continue;}
            else if (appartientppi(resname,"HIS, HID, HIE, PHE, TRP, TYR")){ atm.props.setAromatic(true);continue;}}
        if ((name.compare("CG1")==0 && appartientppi(resname,"ILE, VAL"))
                || (name.compare("CG2")==0 && appartientppi(resname,"ILE, VAL"))){ atm.props.setHydrophobic(true);continue;}
        if (name.compare("CH2")==0 && resname.compare("TRP")==0) {atm.props.setAromatic(true);continue;}
        if (name.compare("CZ")==0){
            if (appartientppi(resname,"PHE, TYR"))                    { atm.props.setAromatic(true);continue;}}
        if (appartientppi(name, "CZ2, CZ3") && resname.compare("TRP")==0)          { atm.props.setAromatic(true);continue;}
        if (name.compare("N")==0)                               { atm.props.setDonor(true);continue;}
        if (appartientppi(name,"ND1, ND2") && appartientppi(resname,"HIS, HID, HIE"))          { atm.props.setDonor(true);atm.props.setAcceptor(true);continue;}
        if (name.compare("ND2")==0 && resname.compare("ASN")==0)                { atm.props.setDonor(true);continue;}
        if (name.compare("CZ")==0 && appartientppi(resname,"ARG"))            { atm.props.setCation(true);continue;} //% Change NE en CZ
        if (name.compare("NE1")==0 && resname.compare("TRP")==0)                { atm.props.setDonor(true);continue;}
        if (name.compare("NE2")==0  ){
            if (resname.compare("GLN")==0)                          {atm.props.setDonor(true);continue;}
            else if (appartientppi(resname,"HIS, HID, HIE"))                   { atm.props.setDonor(true);
                atm.props.setAcceptor(true);continue;}}
        if (name.compare("NZ")==0 && resname.compare("LYS")==0)            {atm.props.setCation(true);continue;}
        if (name.compare("O")==0 && resname.compare("HOH")==0)                               {
            atm.props.setAcceptor(true);atm.props.setDonor(true);continue;}
        else if (name.compare("O")==0){
            atm.props.setAcceptor(true);continue;}

        if (appartientppi(name,"OD1, OD2") && resname.compare("ASP")==0)           {atm.props.setAnion(true);continue;}
        if (name.compare("OD1")==0 && resname.compare("ASN")==0)                {atm.props.setAcceptor(true);continue;}
        if (name.compare("OE1")==0 && resname.compare("GLN")==0)                {atm.props.setAcceptor(true);continue;}
        if (appartientppi(name,"OE1, OE2")  && resname.compare("GLU")==0)          {atm.props.setAnion(true);continue;}
        if ((name.compare("OG") ==0 && resname.compare("SER")==0)||
                (name.compare("OG1")==0 && resname.compare("THR")==0)||
                (name.compare("OH") ==0 && resname.compare("TYR")==0))              {atm.props.setDonor(true);
            atm.props.setAcceptor(true);continue;}
        if (name.compare("OXT")==0)                             {atm.props.setAnion(true);continue;}
        if (appartientppi(name, "SD, SG")){
            if (appartientppi(resname, "CYS, CYX")) {atm.props.setAcceptor(true);continue;}
            else if (resname.compare("MET")==0) {atm.props.setHydrophobic(true);continue;}}
        if (appartientppi(resname, "FE, MG, MG, MN, ZN, CO, CA, NA, ")) {atm.props.setCation(true),atm.props.setMetal(true);continue;}
        atm.props = p;

    }

}


void PPI::detectPPi(std::vector<std::string> &vecChain,int nchaint,bool outputall) // Process
{

    // Interface Proteine Proteine

    //    w_chain = w_chaina;
    //    w2_chain= w_chainb;
    //    ChainSele =chainsele;
    //    ChainSele2=chainseleb;
    nChaintot = nchaint;
    vecChains = vecChain;
    nchain = vecChains.size();

    Grid grid;
    grid.setWHydrogen(true);

    ////Detecte le nombre de chaines et le nombre d'interfaces maximum possible
    PPI::mapinitialisation();

    //    if (mapinteraction.size()==0 )return;
    ////Determine les atomes proches entre toutes les chaines
    PPI::mapdetection();

    //    //% Comment traiter les interface, toute? la plus grande?...
    //% Qu'est ce que c'est moche
//    Molecule &moletest = *(cp_orig.getMole(MoleType::PROTEIN));
//    PPI::defPropsppi(moletest);



    ////Remplis le vecteur interface avec les données Intergene
    PPI::makeinterface();
    //    PPI::makeinterface(volsite);


    sort(intersize.begin(),intersize.end(),comparaisoninterface);

    PPI::calcenfouiss(outputall);

    PPI::detectcavity(outputall);

#ifdef ICHEM_DEBUG
    cout << PDB_name<<"|DETPPI|CHNTRT\t "<< nChaintot<<" chaines traitées " << endl;
#endif


    mapinteraction.clear();

    //    interfaces.clear();
    return;


}
