#include "headers/ICTools/switch.h"

using namespace std;
using namespace ICMole;

void IChemSwitch::helpBSAcalc() const
{

}




void IChemSwitch::BSAcalc() const 
{
    const size_t InputSize = Input_Values.size();
    if (InputSize != 2)
        throw MoleExcept(9020101,
                         "IChem::BSACalc",
                         "Number of parameters must be 2");

    // USER INPUTS :
    const std::string& fProtein = Input_Values.at(0);
    const std::string& fLigand  = Input_Values.at(1);
    string fname="";
    bool multimol2 = true;
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

        if (!multimol2){
            // STEP 1 - LOADING DATA :
            // 1.1 - Loading protein :
            iread.loadNewFile(fProtein);
            iread.loadInComplex(icomplex,MoleType::PROTEIN);
            if (icomplex.getMole(MoleType::PROTEIN)==(Molecule*)NULL)
                throw MoleExcept(9020102,
                                 "IChem::BSACalc",
                                 "No protein found in "+fProtein);
            // 1.2 - Loading Ligand :
            iread.loadNewFile(fLigand);
            iread.loadInComplex(icomplex,MoleType::LIGAND);
            if (icomplex.getMole(MoleType::LIGAND)==(Molecule*)NULL)
                throw MoleExcept(9020102,
                                 "IChem::BSACalc",
                                 "No ligand found in "+fLigand);



            // 1.3 - Writing simplification :
            Molecule &ligand = *icomplex.getMole(MoleType::LIGAND);
            Molecule &protein =*icomplex.getMole(MoleType::PROTEIN);


            // STEP 2 - GENERATING GRID :
            Grid grid;
            grid.setProjLength(2);
            // Create a molecular grid surrounding the protein.
            // Each cube is 0.5Angstroems edge.
            grid.createMolecularBoxes(ligand,0.5);

            // STEP 3 - FINAL VALUES :
            double volLigand=0, BSA=0;
            if (fname.empty()) fname = ligand.getName();
            grid.calcBSA(protein,ligand,BSA,volLigand);

            cout <<fname<<"\t"<< BSA<<"\t"<<volLigand<<endl;

        }
        else{
            // STEP 1 - LOADING DATA :
            // 1.1 - Loading protein :
            iread.loadNewFile(fProtein);
            iread.loadInComplex(icomplex,MoleType::PROTEIN);
            if (icomplex.getMole(MoleType::PROTEIN)==(Molecule*)NULL)
                throw MoleExcept(9020102,
                                 "IChem::BSACalc",
                                 "No protein found in "+fProtein);


            iread.loadNewFile(fLigand);


            const size_t nLigandR = iread.getNumMolecules();

            Molecule &protein =*icomplex.getMole(MoleType::PROTEIN);

            // STEP 2 - GENERATING GRID :
            Grid grid;
            grid.setProjLength(2);


            //STEP 3 READ LIGANDS
            while (!iread.isEOF())
            {

                Molecule ligand;
                iread.loadNextMolecule(ligand,MoleType::LIGAND);
                ligand.checkMOL2();
                grid.createMolecularBoxes(ligand,0.5);

                // STEP 3 - FINAL VALUES :
                double volLigand=0, BSA=0;
                if (fname.empty()) fname = ligand.getName();
                grid.calcBSA(protein,ligand,BSA,volLigand);
                cout <<fname<<"\t"<< BSA<<"\t"<<volLigand<<endl;
            }
        }
    }
    catch (MoleExcept &e)
    {
        e.addTrace("IChem::BSACalc");
        throw;
    }
}



