#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/fragments.h"
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
        << endl
        << "Examples:" << endl
        << "  IChem utils bsa protein.mol2 ligand.mol2" << endl
        << "  IChem utils frag protein.mol2 ligand.mol2" << endl
        << endl
        << "---------------------------------------------------------------------------" << endl
        << endl;
}


void IChemSwitch::utils() const 
{
    const size_t InputSize = Input_Values.size();
    if (InputSize != 3)
        throw MoleExcept(9020101,
                         "IChem::UTS",
                         "Number of parameters must be 3");

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
    }
    catch (MoleExcept &e)
    {
        e.addTrace("IChem::BSACalc");
        throw;
    }



}
