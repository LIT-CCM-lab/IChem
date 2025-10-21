#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/fragments.h"
#include "headers/ICCalcs/pdbconvert.h"
using namespace std;
using namespace ICMole;

void IChemSwitch::helpUtils() const {
    cout
        << "utils - General utilities" << endl
        << "Usage: IChem utils [subcommands] protein.mol2 ligand.mol2" << endl
        << endl
        << "Subcommands:" << endl
        << "  bsa    Compute buried surface area between protein and ligand" << endl
        << "  frag   Fragment ligand within protein context" << endl
        << "  water  Extract water molecules from a PDB file" << endl
        << endl
        << "Examples:" << endl
        << "  IChem utils bsa protein.mol2 ligand.mol2" << endl
        << "  IChem utils frag protein.mol2 ligand.mol2" << endl
        << "  IChem utils water protein.mol2 ligand.mol2 protein.pdb" << endl
        << endl
        << "---------------------------------------------------------------------------" << endl
        << endl;
}


void IChemSwitch::utils() const throw(MoleExcept)
{
    const size_t InputSize = Input_Values.size();
    if (InputSize < 3 || InputSize > 4)
        throw MoleExcept(9020101,
                         "IChem::UTS",
                         "Number of parameters must be either 3 or 4");

    // USER INPUTS :
    const std::string& fJob = Input_Values.at(0);
    const std::string& fProtein = Input_Values.at(1);
    const std::string& fLigand  = Input_Values.at(2);
    string fname="";
    bool multimol2 = false;

    if (!Opt_Values.empty())
        for ( std::map<std::string,std::vector<std::string> >::const_iterator ItOption= Opt_Values.begin();
              ItOption!=Opt_Values.end();
              ++ItOption)
        {
            const std::string& key=(*ItOption).first;

            const std::vector<std::string>& val=(*ItOption).second;

            if (key.compare("-name")==0) fname=val.at(0);
            else if (key.compare("--multim2")==0) multimol2 = true;

        }

    try{
        // PARAMETERS NEEDED FOR OUR WORK :
        Complex icomplex;
        MoleReader iread;


        // STEP 1 - LOADING DATA :
        // 1.1 - Loading protein :
        iread.loadNewFile(fProtein);
        iread.loadInComplex(icomplex,MoleType::PROTEIN);
        if (icomplex.getMole(MoleType::PROTEIN)==(Molecule*)NULL)
            throw MoleExcept(9020102,
                             "IChem::BSACalc",
                             "No protein found in "+fProtein);
        if (!multimol2){
            // 1.2 - Loading Ligand :
            iread.loadNewFile(fLigand);
            iread.loadInComplex(icomplex,MoleType::LIGAND);
            if (icomplex.getMole(MoleType::LIGAND)==(Molecule*)NULL)
                throw MoleExcept(9020102,
                                 "IChem::BSACalc",
                                 "No ligand found in "+fLigand);



            // 1.3 - Writing simplification :
            //            Molecule &ligand = *icomplex.getMole(MoleType::LIGAND);
        }
        Molecule &protein =*icomplex.getMole(MoleType::PROTEIN);

        if (fJob == "bsa")
        {
            if (!multimol2){
                Molecule &ligand = *icomplex.getMole(MoleType::LIGAND);

                // STEP 2 - GENERATING GRID :
                Grid grid;
                grid.setProjLength(2);
                // Create a molecular grid surrounding the protein.
                // Each cube is 0.5Angstroems edge.

                grid.createMolecularBoxes(ligand,0.5);

                // STEP 3 - FINAL VALUES :
                double volLigand=0, BSA=0;
                if (fname.empty()) fname = protein.getName();
                grid.calcBSA(protein,ligand,BSA,volLigand);

                cout <<fname<<"\t"<< BSA<<"\t"<<volLigand<<endl;
            }
            else {
                iread.loadNewFile(fLigand);


                const size_t nLigandR = iread.getNumMolecules();

                Molecule &protein =*icomplex.getMole(MoleType::PROTEIN);

                // STEP 2 - GENERATING GRID :
                Grid grid;
                grid.setProjLength(2);


                //STEP 3 READ LIGANDS
                                    int nb = 0;
                while (!iread.isEOF())
                {

                    Molecule ligand;
                    iread.loadNextMolecule(ligand,MoleType::LIGAND);
                    ligand.checkMOL2();
                    if (nb ==0) {
                        grid.createMolecularBoxes(ligand,0.5);
                    }
                    nb ++;
                    // STEP 3 - FINAL VALUES :
                    double volLigand=0, BSA=0;
                    fname = ligand.getName();
                    grid.calcBSA(protein,ligand,BSA,volLigand);
                    cout <<fname<<"\t"<< BSA<<"\t"<<volLigand<<endl;
                }
            }
        }
        else if (fJob== "frag")
        {
            FragGenerator FragGen(icomplex);
            FragGen.fragmentLigand();
            FragGen.calcInteractions();
            if (fname.empty()) fname = protein.getName();
            FragGen.saveFragments(fname);
        }
        else if (fJob == "water")
        {
            Molecule &ligand = *icomplex.getMole(MoleType::LIGAND);

            const std::string& fPDB  = Input_Values.at(3);
            PDBConvert PDConv(protein.getName());

            PDConv.loadPDB(fPDB);
            Complex pdbcplx;

            try
            {

                PDConv.harmonizeSizeLine();
                PDConv.changeMSEtoMET();
                PDConv.changeCSEtoCYS();
                PDConv.moveHETATMtoend();
                PDConv.selAltAtm();
                PDConv.renumAtms();
                PDConv.toComplex(pdbcplx);

            }
            catch (MoleExcept &e)
            {
                cerr << e.getCode()<<endl
                     << e.getData()<<endl
                     << e.getTrace()<<endl;
                return;
            }


            try
            {



                Molecule::loadRules();
                Residu::loadRules();
                for (unsigned int I=0;I < NB_MOLETYPE;I++)
                {
                    for (unsigned int J=0;J<NB_RESTYPE;J++)
                    {
                        Residu::Rules[I][J]=MoleType::PROTEIN;
                    }
                }

                if (pdbcplx.getMole(MoleType::PROTEIN)==(Molecule*)NULL) throw MoleExcept(999999,"No protein found","");

                Molecule &protein=*pdbcplx.getMole(MoleType::PROTEIN);

                protein.setResiduTypes();

                //protein.selChains(selChain,PDB_ID+"|SELCHA|");

                try{
                    pdbcplx.genGrid(1.5,true);
                }catch (MoleExcept &e)
                {
                    pdbcplx.clearGrid();
                }



                AtomList failedAtom;
                ResiduList failedResidu;
                protein.createConnect(failedResidu, failedAtom,true);
                protein.cleanUnwanted();

                protein.matchTemplate(failedResidu, failedAtom);

                protein.addHydrogen();



                pdbcplx.clearGrid();

                pdbcplx.splitMoleculeInToComplex(protein);

                pdbcplx.selectWater() ;
                double bestdist=1000;
                MoleList tomove;
                for (ItCMole iMole= pdbcplx.firstMole(MoleType::WATER);iMole != pdbcplx.lastMole(MoleType::WATER);++iMole)
                {
                    Molecule& mole=**iMole;
                    bestdist=1000;
                    for (size_t iWatAtm=0;iWatAtm != mole.numAtom();++iWatAtm)
                    {
                        const Atom& atmWat=mole.getAtom(iWatAtm);
                        if (!atmWat.isHydrogen())continue;
                        for (size_t iLigAtm=0;iLigAtm != ligand.numAtom();++iLigAtm)
                        {
                            const Atom& atmLig=ligand.getAtom(iLigAtm);
                            if (atmLig.isHydrogen())continue;
                            const double dist=atmLig.fixpos.calcDist(atmWat.fixpos);
                            if (dist < bestdist)bestdist=dist;
                        }
                        if (bestdist < 6.5)break;
                    }
                    if (bestdist >=6.5)continue;
                    tomove.push_back(&mole);
                }

                if (tomove.empty()) {cerr << "NO WATER FOUND NEAR GIVEN LIGAND"<<endl;return;}
                for (ItMole it = tomove.begin(); it != tomove.end();++it)
                {
                    Molecule& mole = **it;
                    pdbcplx.moveMoleToComplex(&mole,icomplex);

                }

                Residu::loadRules(true);
                ligand.setUse(false);
                MoleWriter mw("output.mol2");mw.setOnlyUsed(true);mw.writeMOL2(icomplex,icomplex.getMole(MoleType::PROTEIN)->getName());
            }catch  (MoleExcept &e)
            {
                cerr << e.getCode()<<endl
                     << e.getData()<<endl
                     << e.getTrace()<<endl;
                return;
            }

        }
    }
    catch (MoleExcept &e)
    {
        e.addTrace("IChem::BSACalc");
        throw;
    }



}
