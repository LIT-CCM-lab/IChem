#include "headers/ICTools/switch.h"

using namespace std;
using namespace ICMole;


void IChemSwitch::helpRealign() const
{
  cout << "realign - Molecular alignment"                                        <<endl
  <<"realign rigidM mobilM applied1 applied2 ....                           "<<endl
  <<" ||-> rigidM : reference molecule to apply alignment to                 "<<endl
  <<" ||-> mobilM : comparison molecule to apply alignment from              "<<endl
  <<" ||-> applied: molecule to apply rotation/translation to                "<<endl
  << endl
  << "  [General options]"<<endl
  << "      -gmatch N (NAME)  Use graph matching to align"                        <<endl
  << "         NAME           Atom Name matching                                 "<<endl
  << "         ATMN           Atomic Name matching                               "<<endl
  << "         MOL2           MOL2 Type matching                                 "<<endl
  << "         CALP           CAlpha Atom matching (protein only)                "<<endl
  << "      --wMob            Also output the aligned mobilM                     "<<endl
  << "      -rule   R         "<<endl
  << " By default, the program will perform an atom by atom match, without taking"<<endl
  << " care of what kind of atom it match. If you want to perform a match by "    <<endl
  << " regarding only some atoms, this index_string is here to do so "            <<endl
  << "       ex : -i '2-3|1-6|23-160'"                                            <<endl
  << "       Will match the second atom from the reference with the third from "  <<endl
  << "       the comparison, the first with the sixth ..."                        <<endl
  <<endl
  << "###########################################################################"<<endl
  <<endl;
}

void IChemSwitch::realign() const throw(MoleExcept)
{
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize < 2) throw MoleExcept(9010601,"IChem::Rotamole","Not enough parameters");

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// PARAMETERS DEFINITION /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Command line Parameters
    const std::string p_rigidF = Input_Values.at(0);  // Rigid molecular file path
    const std::string p_mobilF = Input_Values.at(1);  // Mobile molecular file path
    vector<string> p_toapplyF;                        // List of molecular file path to apply the rotation
    bool p_wMob=false;                                // FALSE: Not outputing mobile aligned molecule
    bool p_localpath=true;
    unsigned int p_gMatch=0;                              // TRUE : use graph matching to align
    std::string p_rule = "";                          // Rule to make alignment
    std::string p_OutPath = "./";                     // Output directory
bool verbose=false;
    // Function parameters
    Complex f_cp;                               // Contains all the molecules to apply alignment
    Complex f_rm;                               // Contains rigid and mobile molecules
    MoleReader f_mread;
    MoleWriter f_mwrite;
    RigidBody f_rb;
    CoordList refCooList,compCooList;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// PARAMETERS ASSIGNATION /////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    for (unsigned int i=2; i<InputSize;i++) p_toapplyF.push_back(Input_Values.at(i));

    {
      string opt_name,opt_val;
      for (std::map<std::string,vector<std::string> >::const_iterator it = Opt_Values.begin(); it != Opt_Values.end(); it++)
        {
          const std::string &opt_name = (*it).first;
          const vector<string> & opt_val = (*it).second;
          const std::string& value = opt_val.at(0);
          //-path --gmatch --wMob -rule
          if (opt_name.compare("-path")     == 0){ p_OutPath=value;p_localpath=false;}
          else if (opt_name.compare("-gmatch")  == 0)
            {
              if (value == "NAME") p_gMatch = 1;
              else if (value == "ATMN") p_gMatch = 2;
              else if (value == "MOL2") p_gMatch = 3;
              else if (value == "CALP") p_gMatch = 4;
              else throw MoleExcept(9010602,"IChem::Rotamole","Unrecognized command for -gmatch option : "+value);
            }
          else if (opt_name.compare("--wMob")    == 0)    p_wMob=true;
          else if (opt_name.compare("-rule")     == 0)
            {
              if (value.find(",")==string::npos || value.find("-")==string::npos)
                throw MoleExcept(9010603,"IChem::Rotamole","Unrecognized rule type : "+value);
              p_rule=value;
            }
          else throw MoleExcept(9010604,"IChem::Rotamole","Unrecognized command : "+opt_name);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// STATUS /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    if (p_gMatch && p_rule.length())throw MoleExcept(9010605,"IChem::Rotamole","--gMatch and -rule cannot be set together");
    if (p_OutPath.substr(p_OutPath.length()-1) != "/") p_OutPath+="/";


    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////// PARAMETERS LOADING ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    if (verbose) cout << "# LOADING FILES"<<endl;;

    /////////////////
    //// REF/MOB ////
    /////////////////
    // Loading molecules :
    f_mread.loadNewFile(p_rigidF); f_mread.loadInComplex(f_rm,MoleType::OTHER);
    f_mread.loadNewFile(p_mobilF); f_mread.loadInComplex(f_rm,MoleType::OTHER);
    // Getting address :
    Molecule &f_rigidM = f_rm.getMoleAtPos(0);
    Molecule &f_mobilM = f_rm.getMoleAtPos(1);
    f_mwrite.setWithRotPos(true);

    /////////////////
    //// APPLIED ////
    /////////////////
    for(vector<string>::iterator it = p_toapplyF.begin(); it != p_toapplyF.end(); it++)
      {
        try
        {
          f_mread.loadNewFile(*it);
          f_mread.loadInComplex(f_cp,MoleType::OTHER);
        }
        catch (MoleExcept &e)
        {
          cerr <<" Unable to load "<< (*it)<<" file "<<endl
              <<" Cause : "<< e.getCode()<<"\t"<<e.getSource()<<"\t"<<e.getData()<<endl;

        }
      }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// CALCULATION ROTATION //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (f_rigidM.numAtom()<3) throw MoleExcept(9010606,"IChem::Rotamole","Rigid molecule must have at least 3 atoms to perform alignement");
    if (f_mobilM.numAtom()<3) throw MoleExcept(9010607,"IChem::Rotamole","Mobile molecule must have at least 3 atoms to perform alignement");

    if (!p_gMatch && p_rule.length()==0)
      {
        if (verbose) cout << "# ALIGNING MOLECULES "<<endl;
        if (f_mobilM.numAtom()<f_rigidM.numAtom()) throw MoleExcept(9010608,"IChem::Rotamole","Mobile molecule must have the same number of atoms as the reference");
        refCooList.reserve(f_rigidM.numAtom());
        compCooList.reserve(f_mobilM.numAtom());
        for (size_t iAtmR=0;iAtmR< f_rigidM.numAtom();++iAtmR)
        {

            refCooList.push_back(f_rigidM.getAtom(iAtmR).fixpos);
        }
        for (size_t iAtmC=0;iAtmC< f_mobilM.numAtom();++iAtmC)
        {
            compCooList.push_back(f_mobilM.getAtom(iAtmC).fixpos);
        }
        f_rb.loadCoordsToRigid(refCooList);
        f_rb.loadCoordsToMobile(compCooList);

      }
    else if (p_gMatch)
      {
        if (verbose) cout << "# COMPUTING GRAPH MATCHING "<<endl;
        CliqueData cliqueData;

        bool rule=false;
        double dist;

        CoordList rlist,clist;
        for (size_t iAtmR=0;iAtmR< f_rigidM.numAtom();++iAtmR)
        {
            const Atom& atmR = f_rigidM.getAtom(iAtmR);
            for (size_t iAtmC=0;iAtmC< f_mobilM.numAtom();++iAtmC)
            {
                const Atom& atmC = f_mobilM.getAtom(iAtmC);
              rule=false;
              switch (p_gMatch)
                {
                case 1: if (atmR.getName()      .compare(atmC.getName())      ==0) rule=true;break;
                case 2: if (atmR.getAtomicName().compare(atmC.getAtomicName())==0) rule=true;break;
                case 3: if (atmR.getMOL2Type()  .compare(atmC.getMOL2Type())  ==0) rule=true;break;
                case 4: if (atmR.getName()=="CA" && atmC.getName()=="CA"
                         &&atmR.getResiduName()==atmC.getResiduName())          rule=true;break;
                }
              if (!rule)continue;
              Pair new_pr(&atmR,&atmC,&cliqueData.graph.addVertex(),1);
//cout << atmR.getIdentifier()<<"\t"<<atmC.getIdentifier()<<endl;
//              new_pr.Id=listPairs.size();
              cliqueData.listPairs.push_back(new_pr);
            }
        }

        //if (verbose)
            cout << "# NUMBER OF POSSIBLE MATCHES : "<< cliqueData.listPairs.size()<<endl;
        for (vector<Pair>::const_iterator itR = cliqueData.listPairs.begin();itR != cliqueData.listPairs.end(); itR++){
            const Pair& PeR = *itR;
            for (vector<Pair>::const_iterator itC = itR+1;itC != cliqueData.listPairs.end(); itC++)
              {
                const Pair& PeC = *itC;
                if (PeR.ref==PeC.ref || PeR.comp==PeC.comp)continue;
                dist = fabs(PeR.ref->fixpos.calcDist(PeC.ref->fixpos)-PeR.comp->fixpos.calcDist(PeC.comp->fixpos));

                if (dist < 0.5) cliqueData.graph.addEdge(*PeR.vertex,*PeC.vertex);

              }
          }


        if (verbose) cout << "# DETECTING CLIQUES ";
        GraphMatch GrM(cliqueData);      if (verbose) cout << ".";
        GrM.createMatrix();             if (verbose) cout << ".";

        //cout<< GrM.printMatrix()<<endl;
        GrM.runSearchAllMaxCliques(verbose);   if (verbose) cout << ".";
        if (verbose) cout << " "<<cliqueData.cliques.size() << " cliques found"<<endl;
        vector<int> &bestclique= *cliqueData.cliques.begin();
        for (vector<int>::const_iterator it=bestclique.begin(); it != bestclique.end(); it++)
          {
            if (verbose) cout << " # MATCH : "<< cliqueData.listPairs.at(*it).ref->getIdentifier() <<"\t"<< cliqueData.listPairs.at(*it).comp->getIdentifier()<<endl;
            rlist.push_back(cliqueData.listPairs.at(*it).ref->fixpos);
            clist.push_back(cliqueData.listPairs.at(*it).comp->fixpos);
          }
        f_rb.loadCoordsToRigid(rlist);
        f_rb.loadCoordsToMobile(clist);

      }
    else if (p_rule.length())
      {
        vector<string> tokens;
        unsigned int RI=0,CI=0;
        CoordList rlist,clist;
        tokenStr(p_rule,tokens,",");
        for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++)
          {
            RI=atoi((*it).substr(0,(*it).find("-")).c_str())-1;
            CI=atoi((*it).substr((*it).find("-")+1).c_str())-1;
            if (verbose) cout << "# "<<RI<<":"<<f_rigidM.getAtom(RI).getIdentifier()<<" with "<< CI<<":"<<f_mobilM.getAtom(CI).getIdentifier()<<endl;
            if (RI > f_rigidM.numAtom()) throw MoleExcept(9010609,"IChem:Rotamole","Given atom id is above the number of atom in the reference molecule for rule :"+(*it));
            if (CI > f_mobilM.numAtom()) throw MoleExcept(9010610,"IChem:Rotamole","Given atom id is above the number of atom in the comparison molecule for rule :"+(*it));
            rlist.push_back(f_rigidM.getAtom(RI).fixpos);
            clist.push_back(f_mobilM.getAtom(CI).fixpos);
          }

        f_rb.loadCoordsToRigid(rlist);
        f_rb.loadCoordsToMobile(clist);

      }

    double rmsd = f_rb.calcRotation();
    if (verbose) cout << "# FINAL RMSD : " <<rmsd<<endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// MAKING ALIGNMENT ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
  unsigned int numM=0;
  string mpath,fname;
  if (p_localpath) mpath+="rot_";
  else mpath+=p_OutPath;

  if (f_cp.getNumMolecule()==0)
    {
  //f_rb.mobilToRef(f_mobilM);
  fname=p_mobilF;
  if (fname.find_last_of("/") != string::npos) fname=mpath+fname.substr(fname.find_last_of("/")+1);
  else fname = mpath+fname;
  if(verbose)cout << "# SAVING "<< fname<<endl;
  f_mwrite.newFile(fname);
  f_mwrite.writeMOL2(&f_mobilM);
    }

    for (size_t iMole=0; iMole < f_cp.getNumMolecule();++iMole)
      {
        Molecule &mole = f_cp.getMoleAtPos(iMole);
        CoordList coords;coords.clear();
        for (size_t iAtm =0; iAtm < mole.numAtom();++iAtm) coords.push_back(mole.getAtom(iAtm).fixpos);
        f_rb.mobilToRef(coords);
        size_t iAtm=0;
        for (ItCAtom itA= mole.firstAtom();itA != mole.lastAtom();++itA)
        {
            (*itA)->rotpos=coords.at(iAtm);
            iAtm++;
        }



        fname=p_toapplyF.at(numM);
        if (fname.find_last_of("/") != string::npos) fname=mpath+fname.substr(fname.find_last_of("/")+1);
        else fname = mpath+fname;
        if(verbose)cout << "# SAVING "<< fname<<endl;
        f_mwrite.newFile(fname);
        f_mwrite.writeMOL2(&mole);
        numM++;
      }
    if (p_wMob)
      {
      //  f_rb.mobilToRef(f_mobilM);
        CoordList coords;coords.clear();
        for (size_t iAtm =0; iAtm < f_mobilM.numAtom();++iAtm) coords.push_back(f_mobilM.getAtom(iAtm).fixpos);
        f_rb.mobilToRef(coords);
        size_t iAtm=0;
        for (ItCAtom itA= f_mobilM.firstAtom();itA != f_mobilM.lastAtom();++itA)
        {
            (*itA)->rotpos=coords.at(iAtm);
            iAtm++;
        }
        fname=p_mobilF;

        if (fname.find_last_of("/") != string::npos) fname=mpath+fname.substr(fname.find_last_of("/")+1);
        else fname = mpath+fname;
        if(verbose)cout << "# SAVING "<< fname<<endl;
        f_mwrite.newFile(fname);
        f_mwrite.writeMOL2(&f_mobilM);
      }


}
