
#include "headers/ICTools/switch.h"
#include "headers/ICMole/similarity.h"
#include "headers/ICCalcs/interaction.h"



using namespace std;
using namespace ICMole;

void IChemSwitch::helpFGPS() const
{
    cout<< "Fingerprint Similarity"                                                     <<endl
        << "    sims ref comp                                                       (1)"<<endl
        << "    sims file                                                           (2)"<<endl
        << "    sims RefInt CompInt                                                 (3)"<<endl
        << "    sims RefInt RefInt CompInt CompInt                                               (4)"<<endl
        << endl
        << "  [General options]"<<endl
        << "      --wInts           To use interactions instead of fingerprints "<<endl
        << "                        For (3) only"<<endl
        << "      --small           Use small fingerprint. For (3) only"<<endl
        << "      --binary          To add when the fingerprint is binary"              <<endl
        << "      -metric N (TC)    Select the metric"                                  <<endl
        << "         TC             Tanimoto metric"                                    <<endl
        << "         HM             Hamming distance"                                   <<endl
        << "         RT             Ref Tversky "                                       <<endl
        << "         FT             Fit Tversky "                                       <<endl
        << "         DI             Dice "                                              <<endl
        << "         SO             Soergel "                                           <<endl
        <<endl
       << "###########################################################################"<<endl;

}


Fingerprint Fingerprint::generateIFP(const std::string& protein_file, const std::string& ligand_file, bool numeric) {
    
    using namespace ICMole;

    Complex complex;
    MoleReader reader;

    reader.loadNewFile(protein_file);
    reader.get_format_file();
    reader.loadInComplex(complex, MoleType::PROTEIN);

    reader.loadNewFile(ligand_file);
    reader.loadInComplex(complex, MoleType::LIGAND);

    Molecule* ligand = complex.getMole(MoleType::LIGAND);
    if (!ligand) {
        throw MoleExcept(9020102, "generateIFPfromFiles", "No ligand found in " + ligand_file);
    }

    ligand->checkMOL2();
    ligand->ringPerception();

    Interactions ints(complex);
    InterResults res;
    ints.detectInteractions(*ligand, res, true, true);  

    ints.genIFP(res, 0);

    res.IFP.setName(ligand->getName());

    return res.IFP;
}


void IChemSwitch::runFGPS() const throw(ICMole::MoleExcept)
{

    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize != 1 && InputSize != 2 && InputSize != 4) {throw MoleExcept(9010401,"IChem::runFGPS","Not enough parameters");}

    bool numeric=true;
    unsigned int metric = 1;
    bool wInts=false;
    bool full=true;


    for (std::map<std::string,std::vector<std::string> >::const_iterator it = Opt_Values.begin(); it != Opt_Values.end(); it++) {
        
        const std::string &opt_name = (*it).first;
        const vector<string> & opt_val = (*it).second;
        const std::string& value = opt_val.at(0);


        if (opt_name.compare("--binary")==0) {numeric=false;}
        else if (opt_name.compare("--wInts")==0){wInts=true;}
        else if (opt_name.compare("--small")==0){full=false;}

        else if (opt_name.compare("-metric")    == 0) {
            if (value == "TC") metric=1;
            else if (value == "HM") metric=2;
            else if (value == "RT") metric=3;
            else if (value == "FT") metric=4;
            else if (value == "DI") metric=5;
            else if (value == "SO") metric=6;
            else throw MoleExcept(9010402,"IChem::runFGPS","Unrecognized option for metric : "+value);
        }


        else throw MoleExcept(9010404,"IChem::runFGPS","Unrecognized optionp : "+opt_name);
    }

    if (InputSize == 2) {

        const std::string Ref = Input_Values.at(0);
        const std::string Comp = Input_Values.at(1);


        string ligne;
        size_t pos;
        vector<Fingerprint> refFGPS;
        vector<Fingerprint> compFGPS;
        Similarity sims(numeric);double simil=0;

        ifstream ifs(Ref.c_str());
        while(!ifs.eof())
        {
            std::getline(ifs,ligne);
            if (ligne.length()==0)continue;
            pos = ligne.find_first_of("\t");

//            std::cout << ligne.substr(pos+1) << std::endl;

            Fingerprint fgp(ligne.substr(pos+1),ligne.substr(0,pos),numeric);

//            std::cout << fgp.toString() << std::endl;

            refFGPS.push_back(fgp);
        }

        ifs.close();

        ifstream ifsc(Comp.c_str());
        while(!ifsc.eof())
        {
            std::getline(ifsc,ligne);
            if (ligne.length()==0)continue;
            pos = ligne.find_first_of("\t");

            Fingerprint fgpc(ligne.substr(pos+1),ligne.substr(0,pos),numeric);
            compFGPS.push_back(fgpc);
        }

        ifsc.close();

        for (vector<Fingerprint>::iterator it = refFGPS.begin(); it != refFGPS.end(); it++) {

            sims.setRef(*it);
            for (vector<Fingerprint>::iterator itC = compFGPS.begin(); itC != compFGPS.end(); itC++) {
                sims.setComp(*itC);

                switch(metric)
                {
                case 1:simil=sims.Tanimoto();break;
                case 2:simil=sims.Hamming(); break;
                case 3:simil=sims.RTve();    break;
                case 4:simil=sims.FTve();    break;
                case 5:simil=sims.Dice();    break;
                case 6:simil=sims.Soergel(); break;
                }
                cout << (*it).getName()<<"\t"<<(*itC).getName()<<"\t"<<simil<<endl;
            }

        }

    }

    else {
        const std::string FName = Input_Values.at(0);

        string ligne;
        size_t pos;
        vector<Fingerprint> FGPS;
        Similarity sims(numeric);double simil=0;


        ifstream ifs(FName.c_str());
        while(!ifs.eof())
        {
            std::getline(ifs,ligne);
            if (ligne.length()==0)continue;
            pos= ligne.find_first_of("\t");

            Fingerprint fgp(ligne.substr(pos+1),ligne.substr(0,pos),numeric);

            FGPS.push_back(fgp);
        }
        ifs.close();


        for (vector<Fingerprint>::iterator it = FGPS.begin(); it != FGPS.end(); it++)
        {
            sims.setRef(*it);
            for (vector<Fingerprint>::iterator itC = FGPS.begin(); itC != FGPS.end(); itC++)
            {
                sims.setComp(*itC);

                switch(metric)
                {
                case 1:simil=sims.Tanimoto();break;
                case 2:simil=sims.Hamming(); break;
                case 3:simil=sims.RTve();    break;
                case 4:simil=sims.FTve();    break;
                case 5:simil=sims.Dice();    break;
                case 6:simil=sims.Soergel(); break;
                }
                cout << (*it).getName()<<"\t"<<(*itC).getName()<<"\t"<<simil<<endl;
            }

        }

    }

}

