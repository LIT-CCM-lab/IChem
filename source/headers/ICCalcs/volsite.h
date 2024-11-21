#ifndef VOLSITE_H
#define VOLSITE_H
#include <iostream>
//#define ICHEM_DEBUG
#include <math.h>
#include "headers/ICMole/grid.h"
//% J'aime pas trop ca moi
#include "headers/ICSVM/svm_common.h"
#include "headers/ICMole/molecule.h"
#define PTS 40000
#define K 50

namespace ICMole
{

struct VolSiteResult {
    Molecule moleRes;

    unsigned short type, NumCav;
    double drugg,ligRecov,cavRecov,volume;
    int CA,CZ,O,OD1,N,NZ,OG,DU, nbcube ;
    VolSiteResult() :drugg(-1),ligRecov(0),cavRecov(0),volume(0),CA(0),CZ(0),O(0),OD1(0),N(0),NZ(0),OG(0),DU(0),nbcube(0) {}

    void ToString(){
        std::cout << "############ VolSite Result :" << moleRes.getName() << "############" << std::endl
                  << "type : " << type << "\t NumCav : " << NumCav << std::endl
                  << moleRes.getName().substr(0,4) << "|CAVITY|DESCAV\t" << volume <<"\t"<< nbcube <<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<std::endl;
    }



};

struct ClusterPha
{
    Coords Centr;
    unsigned int pha;
    //        HBDONOR=0,
    //        HBACCEPTOR=1,
    //        HYDROPHOBIC=2,

    unsigned int idClus;

};



class VolSite
{
    friend class Complex;
#ifndef ICHEM_DEBUG
    Grid &grid;
#endif
    std::vector<unsigned int> projs;
    /*!< \brief Minimal projection value to consider cube as buried. Default : 40 */
    int  minProj,numCav;
    bool proteinOnly;
    double min_size;
    double max_size;
    double drug_value;
    int proj_lens;
    std::string  name;                                                                                    /*!< \brief Name of cavities */
    BoxList  cubes_cavite;                                                                            /*!< \brief All cubes of cavity */
    double L1,L2,L3,L4,L5,L6,L7;
    std::vector<ClusterPha> clusterList;
    std::vector<ClusterPha>::iterator BestCl;
    int cavitykeepP;
    int cavitykeepC;
    bool tolerance;
    bool Bweight;
    double PHA_N;
    double PHA_NZ;
    double PHA_O;
    double PHA_OD1;
    double PHA_CA;
    double PHA_CZ;
    double PHA_Zn;
    double dist_noH;
    double dist_h;

    //     const  double Angl_H = M_PI;        //180
    //     const double AngT_H = M_PI/3;      //40
    //     const double Angl_D = M_PI/1.65;   //109
    //     const double D120         = M_PI/1.5;    //120
    //     const double D180         = M_PI;
    //     const double angl_aro     = M_PI/2;  // 90

    std::multimap<unsigned int,unsigned int> size_cavs;
    std::map<unsigned int, std::set<const Box*>  > cavities;
    std::set<const Box*> pha_cavity;


public:
    bool clusterPHA_succes;

#ifdef ICHEM_DEBUG
    Grid &grid;
#endif
    VolSite(Complex &cp, const double &boxStep=1.5, const int proj_lens=8, bool def_alternative =false) throw(MoleExcept); // ppi
    VolSite(Complex &cp, Molecule &ligand,   const double &boxStep=1.5, const Coords& box_size=Coords(20,20,20), const int proj_lens=8, bool def_alternative =false) throw(MoleExcept);

    VolSite(Molecule &protein, Molecule &ligand, Grid &grid, const double &boxStep=1.5, const Coords& box_size=Coords(20,20,20), const int proj_lens=8, bool def_alternative =false) throw(MoleExcept); //voslite
    VolSite(Molecule &molecule, Grid &grid, const double &boxStep=1.5, const int proj_lens=8, bool def_alternative =false) throw(MoleExcept); // scoring

    void treatments(unsigned int minNeighbors=9,const double dist_prot=2.5,bool siteout=false) throw(MoleExcept);
    void KmeansPHA(int nb_clust, bool prop=false);
    void clusterCAV(const int &min_proj=40);
    void clusterPHA(int i, int j, int k);
    void proj(const int &min, const bool& force_stop=false);
    void cavToMOL(Molecule& mollen,std::vector<const Box *> &boxlist, const std::string &molename, const std::string &size) const;
    void cavToMOLB(Molecule& molen,  std::vector<const Box*> &boxlist, const std::string& molename, const std::string& size,const int bestb) const;
    void cavToPHA(Molecule& mollen, PharmList &MonPharma, std::vector<const Box *> &boxlist, const std::string &molename, const std::string &size, const int mod, int& bestb) const;
    void cubeToPHA(Molecule& mollen, PharmList &MonPharma,Molecule& ligand, const std::string &molename, const std::string &size, const int mod, const bool check=false, const bool exclu=false) const;
    void defProps(Molecule& protein, bool pharma = false);
    void output_cavs_lig(const Molecule &ligand, const std::string& name, const bool &svm, double &dist, bool Pharma, int PHA_mod);
    void output_cavs_all(const std::string& name, const bool &svm,bool Pharma, int PHA_mod, std::vector<std::string> Input);
    void output_cavs_ppi(const std::string& name, const bool &svm);
    void getDruggability(double &dist, Molecule const & mole, const bool &svm=false, int combo=0, int numcav=0);
    void read_VolSitemodel(MODEL &modelin, const std::string msvm);
#ifdef ICHEM_DEBUG
    void printProj(const std::string& fname) const;
#endif

    void structCavLig(VolSiteResult result[], const Molecule &ligand, const std::string& name, const  bool& svm);
    void structCavAll(VolSiteResult result[], const std::string& name, const bool &svm);


    inline void setMinSize(const unsigned int& min) { min_size = min;}
    inline void setMaxSize(const unsigned int& max) { max_size = max;}
    inline void setDruggabilty(const unsigned int& t) {drug_value=t;}
    inline void setCavityKP(const unsigned int& min) { cavitykeepP = min;}
    inline void setCavityKC(const unsigned int& t) {cavitykeepC=t;}
    inline void setTolerance(const bool& t) {tolerance=t;}
    inline void setWeigth(const bool& t) {Bweight=t;}



    inline void setproteinOnly(bool set) {proteinOnly=set;}
    inline void setName(std::string nom) {name=nom;}
    inline bool getPhaSuccess() const { return clusterPHA_succes;}
    inline void setPhaSuccess(bool success) { clusterPHA_succes=success;}

    const int& getNumCav() const { return numCav;}
    const Grid& getGrid() const {return grid;}
    const unsigned int& getProjValue(const Box& bx) const;
};

}

#endif // VOLSITE_H
