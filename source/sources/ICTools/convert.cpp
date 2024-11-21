    #include "headers/ICTools/switch.h"
    #include "headers/ICCalcs/pdbconvert.h"
    using namespace std;
    using namespace ICMole;

    void IChemSwitch::helpPDBConvert() const
    {
        cout
                << "PDB Process                                                                     "<<endl
                << "     pdbconv protein[.pdb|.mol2] output_dir pdb_id                              "<<endl
                << endl
                << "   --wMOL2    : Use MOL2 File as Input. PDB Options are not available           "<<endl
                << "   --wUnDrug  : Output undruggable cavities"                                     <<endl
                << "   --noLig    : PDB with no Ligand"                                              <<endl
                << "By default all the following options are included.                              "<<endl
                << "All chains will be kept                                                         "<<endl
                << "  [PDB Options]                                                                 "<<endl
                << "   --HARMSIZE : Harmonize size line to 80 characters                            "<<endl
                << "   --MSEMET   : Change MSE to MET                                               "<<endl
                << "   --CSECYS   : Change CSE to CYS                                               "<<endl
                << "   --MOVHET   : move HETATM to the end of file                                  "<<endl
                << "   --ALTATM   : select alternative atoms                                        "<<endl
                << "   --NUMATM   : renumerotate atoms                                              "<<endl
                << "   --UPDMAS   : update the MASTER line                                          "<<endl
                << "   --TOMOL2   : convert to a molecular representation (instead of flat file)    "<<endl
                << " if you use one of the option below, you MUST use also --TOMOL2 option          "<<endl
                << " or use --wMOL2 option                                                          "<<endl
                << "  [MOL2 Options]                                                                "<<endl
                << "   --RESTYP   : apply Residu Class (cofactor/STD_AA/MOD_AA/Ligand ...)          "<<endl
                << "   --BONDSE   : recreate bonds                                                  "<<endl
                << "   --CLNUNW   : clean unwanted residus                                          "<<endl
                << "   --MOL2TY   : apply MOL2 types according to templates                         "<<endl
                << "   --SPLITM   : split molecule into protein/ligand/solvent                      "<<endl
                << "   -SelChain N : List of chains to keep, separated by underscore                "<<endl
                << "   --SELWAT   : select water molecules                                          "<<endl
                << "   --SELLIG   : select ligand                                                   "<<endl
                <<endl
               << "################################################################################"<<endl
               <<endl;

    }
    void IChemSwitch::PDBConversion() const    throw(ICMole::MoleExcept)
    {
        const unsigned int InputSize = (unsigned int)Input_Values.size();
        if (InputSize != 3) {throw MoleExcept(9010301,"IChem::pdbconv","Not enough parameters");}

        const string PDBFile = Input_Values.at(0);
        const string Dir = Input_Values.at(1);
        const string PDB_ID = Input_Values.at(2);

        unsigned int NRule=0;
        bool harmsize=false;
        bool msemet=false;
        bool csecys=false;
        bool movhet=false;
        bool altatm=false;
        bool numatm=false;
        bool updmas=false;
        bool tomol2=false;
        bool restyp=false;
        bool bondse=false;
        bool clnunw=false;
        bool mol2ty=false;
        bool splitm=false;
        bool selwat=false;
        bool sellig=false;
        bool full  =false;
        bool wMOL2 =false;
        bool nolig = false;
        bool keepundrug=true;
        bool addH=false;
        bool process =false;
        string selChain="";

        for (std::map<std::string,vector<std::string> >::const_iterator it = Opt_Values.begin(); it != Opt_Values.end(); it++)
        {
            const std::string &opt_name = (*it).first;
            const vector<string> & opt_val = (*it).second;
            const std::string& value = opt_val.at(0);
            //                         PARAM               SETTING FOR REFERENCE COMPLEX                    SETTING FOR COMPARISON COMPLEX
            if (opt_name.compare("--HARMSIZE") == 0) { harmsize=true; NRule++;}
            else if (opt_name.compare("--MSEMET")   == 0) { msemet=true; NRule++;}
            else if (opt_name.compare("--CSECYS")   == 0) { csecys=true; NRule++;}
            else if (opt_name.compare("--MOVHET")   == 0) { movhet=true; NRule++;}
            else if (opt_name.compare("--ALTATM")   == 0) { altatm=true; NRule++;}
            else if (opt_name.compare("--NUMATM")   == 0) { numatm=true; NRule++;}
            else if (opt_name.compare("--UPDMAS")   == 0) { updmas=true; NRule++;}
            else if (opt_name.compare("--TOMOL2")   == 0) { tomol2=true; NRule++;}
            else if (opt_name.compare("--RESTYP")   == 0) { restyp=true; NRule++;}
            else if (opt_name.compare("--BONDSE")   == 0) { bondse=true; NRule++;}
            else if (opt_name.compare("--CLNUNW")   == 0) { clnunw=true; NRule++;}
            else if (opt_name.compare("--MOL2TY")   == 0) { mol2ty=true; NRule++;}
            else if (opt_name.compare("--SPLITM")   == 0) { splitm=true; NRule++;}
            else if (opt_name.compare("--SELWAT")   == 0) { selwat=true; NRule++;}
            else if (opt_name.compare("--SELLIG")   == 0) { sellig=true; NRule++;}
            else if (opt_name.compare("--wUnDrug")  == 0) { keepundrug=true;     }
            else if (opt_name.compare("--wMOL2" )   == 0) {  wMOL2=true;}
            else if (opt_name.compare("-SelChain")   == 0) { selChain = value; }
            else if (opt_name.compare("--process")   == 0) { process = true; }
            else if (opt_name.compare("--addh")   == 0) { addH = true; }
            else if (opt_name.compare("--noLig" )   == 0) {nolig=true;harmsize=true;msemet=true;csecys=true;movhet=true; altatm=true;numatm=true;updmas=true;tomol2=true;restyp=true;bondse=true;clnunw=true;mol2ty=true; selwat=true; NRule++;}
            else throw MoleExcept(9010310,"IChem::runPDB","Unrecognized optionp : "+opt_name);
        }
        if (NRule==0)full=true;

        PDBConvert PDConv(PDB_ID);
        PDBConvert::for_scPDB=true;
        Complex complex;
        if (!wMOL2)
        {

            PDConv.setPDBName(PDB_ID);
            PDConv.loadPDB(PDBFile);

            try
            {

                if (full || harmsize) PDConv.harmonizeSizeLine();
                if (full || msemet  ) PDConv.changeMSEtoMET();
                if (full || csecys  ) PDConv.changeCSEtoCYS();
                if (full || movhet  ) PDConv.moveHETATMtoend();
                if (full || altatm  ) PDConv.selAltAtm();
                if (full || numatm  ) PDConv.renumAtms();
                if (full || tomol2  ) PDConv.toComplex(complex);

            }
            catch (MoleExcept &e)
            {
                cerr << e.getCode()<<endl
                     << e.getData()<<endl
                     << e.getTrace()<<endl;
                return;
            }
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
            if (wMOL2)
            {
                MoleReader mr(PDBFile);
                mr.loadInComplex(complex,MoleType::PROTEIN);
            }
            if (complex.getMole(MoleType::PROTEIN)==(Molecule*)NULL) throw MoleExcept(999999,"No protein found","");

            Molecule &protein=*complex.getMole(MoleType::PROTEIN);

            if (full || restyp) protein.setResiduTypes(PDB_ID+"|RESTYP|");

            if (selChain.length())protein.selChains(selChain,PDB_ID+"|SELCHA|");

            try{
                complex.genGrid(1.5,true);
            }catch (MoleExcept &e)
            {
                complex.clearGrid();
            }



            AtomList failedAtom;
            ResiduList failedResidu;
            if (full || bondse)protein.createConnect(failedResidu, failedAtom,true,PDB_ID+"|CONECT|");
            if (full || clnunw)protein.cleanUnwanted(PDB_ID+"|DELUNW|");


            if (full || mol2ty)protein.matchTemplate(failedResidu, failedAtom,PDB_ID+"|MATCHT|",true);


            if (addH) protein.addHydrogen(PDB_ID+"|ADDHYD|");



            complex.clearGrid();

            if (full || splitm) complex.splitMoleculeInToComplex(protein,PDB_ID+"|SPLITM|");

            if (full && (complex.getNumMolecule(MoleType::LIGAND)+complex.getNumMolecule(MoleType::COFACTOR) == 0))
            {
                MoleWriter mw(PDB_ID+".mol2");
                mw.writeMOL2(complex,PDB_ID);
                throw MoleExcept(99999,"","No ligand found");

            }
            else if ( nolig) {
                MoleWriter mw(PDB_ID+".mol2");
                mw.writeMOL2(complex,PDB_ID);
            }
            if (full || selwat) complex.selectWater(PDB_ID+"|WATERS|") ;

            sortAndUnique(failedResidu);
            sortAndUnique(failedAtom);

            Residu::loadRules(true);
    for (ItCMole it = complex.firstMole(); it != complex.lastMole(); ++it)
        (*it)->setUse(true);
            if (full || sellig)complex.selectLigand(Dir,failedResidu,failedAtom,PDB_ID+"|SELLIG|",keepundrug);

            if (process){
            ofstream ofs;
            ofs.open(".pdbconv",ios::out|ios::app);
            ofs << "Success" <<endl;
            ofs.close();
    }
        }
        catch (MoleExcept &e)
        {
            cerr << e.getCode()<<endl
                 << e.getData()<<endl
                 << e.getTrace()<<endl;
            return;
        }


    }
