#ifndef INTERACTION_H
#define INTERACTION_H


#include "headers/ICMole/molecule.h"
#include "headers/ICMole/fingerprint.h"
#include <config_options.h>


namespace ICMole
{

 struct InterPoint
{
    int point;
    ICMole::Atom* Prot_Ref;
    ICMole::Atom* Lig_Ref;
    Coords center; // center of the two points
    unsigned int interaction;
    double dist;
    double angle;
    int merged_to;

    InterPoint(const InterPoint  &p):
        point(p.point),
        Prot_Ref(p.Prot_Ref),
        Lig_Ref(p.Lig_Ref),
        center(p.center),
        interaction(p.interaction),
        dist(p.dist),
        angle(p.angle),
        merged_to(p.merged_to) {}

    // InterPoint(const  int& p,
    //             Atom* const protein,
    //             Atom* const ligand,
    //            const Coords& center,
    //            const unsigned int& int_type,
    //            const double& dist,
    //            const double& angl= -100000):
    //             point(p),
    //             Prot_Ref(protein),
    //             Lig_Ref(ligand),
    //             center(center),
    //             interaction(int_type),
    //             dist(dist),
    //             angle(angl),
    //             merged_to(-1){}
    InterPoint(const int& p,
                       Atom* const protein,
                       Atom* const ligand,
                       const Coords& center,
                       const unsigned int& int_type,
                       const double& dist,
                       const double& angl = 0)
    : point(p),
      Prot_Ref(protein != nullptr ? protein : nullptr),
      Lig_Ref(ligand != nullptr ? ligand : nullptr),   
      center(center),
      interaction(int_type),
      dist(dist),
      angle(angl),                
      merged_to(-1) {}   


    // friend std::ostream& operator<<(std::ostream& os, const InterPoint& interPoint) {
    //     std::cout << "\n"; 
    //     // os << "Protein_Ref: " << interPoint.Prot_Ref->getName() << std::endl;
    //     // os << "Ligand_Ref: " << interPoint.Prot_Ref->getName() << std::endl;
    //     os << "Distance: " << interPoint.dist << std::endl;
    //     os << "center.x: " << interPoint.center.x << std::endl;
    //     os << "center.y: " << interPoint.center.y << std::endl;
    //     os << "center.z: " << interPoint.center.z << std::endl;
    //     os << "Interaction: " << interPoint.interaction << std::endl;
    //     os << "Angle: " << interPoint.angle << std::endl;
    //     os << "Merged_to: " << interPoint.merged_to << std::endl;
    //     return os;
    // }
};


 struct InterResults
 {
   Molecule Ints;
   unsigned short CA,CZ,O,N,NZ,OD1,Zn;
   std::vector<InterPoint> listInters;
   std::string IFPString;
   Fingerprint IFP;
   std::string chainA;
   std::string chainB;

   
   InterResults():CA(0),CZ(0),O(0),N(0),NZ(0),OD1(0),Zn(0),chainA(""),chainB(""){}

  // friend std::ostream& operator<<(std::ostream& os, InterResults& results) {
  //   bool first = true;
  //   os << "############### BEGINNING ######################### ";
  //   for(const auto& x: results.listInters) {
  //     if(!first) {
  //       os << "------------";
  //     }
  //     os << x;
  //     first = false;
  //   }
  //   std::cout << "\n";
  //   os << "################# END ############################";
  //   return os;
  // }

  // void displayInterResults() {
  //   for(const auto x: listInters) {
  //     cout << "display listInters elem: " << x << endl;
  //   }
  // } 

 };


class Interactions
{
            Complex& complex;
            Grid &grid;
    // Distances modifiables :

            /**
              * @brief HBond maximal distance
              */
             static double  Dist_H;


             /**
              * @brief Hydrophobic maximal distance
              */
             static double  Dist_Hyd;



             /**
              * @brief Ionic maximal distance
              */
             static double  Dist_Ionic;



             /**
              * @brief Metal/Acceptor maximal distance
              */
             static double  Dist_Metal;



             /**
              * @brief Aromatic maximal distance between center
              */
             static double  Dist_Arom;


             /**
              * @brief HBond minimal distance
              */
             static double  dist_H;


             /**
              * @brief Hydrophobic minimal distance
              */
             static double  dist_Hyd;



             /**
              * @brief Ionic minimal distance
              */
             static double  dist_Ionic;



             /**
              * @brief Metal/Acceptor minimal distance
              */
             static double  dist_Metal;



             /**
              * @brief Aromatic minimal distance between center
              */
             static double  dist_Arom;


    // Angles modifiables :

             /**
              * @brief Median angle for HBond interactions
              */
             static double  Angl_H;


             /**
              * @brief +/- possible deviation for HBond angle
              */
             static double  AngT_H;


             /**
              * @brief Median Aromatic Face to Face angle
              */
             static double  Angl_AromFF;


             /**
              * @brief +/- possible deviation for Aromatic Face to Face angle
              */
             static double  AngT_AromFF;


             /**
              * @brief Median Aromatic Edge to Face angle
              */
             static double  Angl_AromEF;

             /**
              * @brief possible deviation for Aromatic Edge to Face angle
              */
             static double  AngT_AromEF;


             /**
              * @brief PI-Cation maximal distance
              */
             static double Dist_PiCation;


             /**
              * @brief Median angle for PI-Cation interaction
              */
             static double Angl_PiCation;


             /**
              * @brief Possible deviation for PI-Cation interaction
              */
             static double AngT_PiCation;


             /**
              * @brief Weak H-Bond maximal distance
              */
             static double Dist_WHBond;


             /**
              * @brief Median angle for Weak H-Bond interaction
              */
             static double Angl_WHBond;


             /**
              * @brief Possible deviation for Weak H-Bond interaction
              */
             static double AngT_WHBond;


             /**
              * @brief Filtering rules to detect only some interactions
              * Set to true when the interaction needed to be detected
              */
             static bool    wInterType[NB_INTTYPE];

             

             /**
             * @brief Set default values for interaction rules:
             */
            static void setStdRules();

            /**
             * @brief Merge
             * hydrophobic interactions
             * @param interResult : interaction set to be merged
             */
            void mergeInteractions(InterResults& interResult)const;

            /**
             * @brief Merge All interactions
             * @param interResult : interaction set to be merged
             */

            void mergeSpeInts(InterResults& interResult)const;

public:
            /**
              * @brief Constructor
              * @param complex with a least one molecule set a protein
              */
             Interactions(ICMole::Complex &complex) throw(ICMole::MoleExcept);


             /**
              * @brief Determine interaction molecule buriedness
              */
             void calcenfouiss(std::string mol_name);




             /**
              * @brief Detect interaction according to geometric rules
              * @param ligand : molecule to detect interaction with
              * @param interResult : Detected interactions are score here
              * @param wMerge : Merge hydrophobic interactions
              */
             void calcInteractions(Molecule&     ligand,
                                  InterResults& interResult,
                                  bool          wMerge=true,
                                  bool          oldh=true,
                                  bool mono_prop=false, 
                                  bool displayProperties = false,
                                  bool out_lig=false,
                                  bool stdout=false
                                  ) const; //bool displayProperties = false

                
            void displayAtomProperties(ICMole::Molecule& ligand, Atom& atomL, std::string molecule_name) const;


             /**
              * @brief Detect interaction according to geometric rules
              * @param interResult : Detected interactions are score here
              * @param wMerge : Merge hydrophobic interactions
              * @param wHydrogen : Check hydrigen interaction
              */
             void calcInteractionsppi(InterResults& interResult,
                                   bool          wMerge=true,
                                   bool wHydrogen=false,
                                   bool          oldh=false) const;


             /**
              * @brief Convert interResults in a human readable string
              * @param interResult : interactions set
              * @return string describing the binding mode
              */
             std::string toString(const  InterResults& interResult)const;


             /**
              * @brief Convert a binding mode into a molecule
              * @param interResults : interactions set to be converted
              * @param Ligand : True when generating ligand positionned molecule
              * @param Protein : True when generating protein positionned molecule
              * @param Center : True when generating centered positionned molecule
              * @param Merge : True when generation all positions molecule
              */
             void interToMOL2(  InterResults& interResult,
                               const bool&Ligand,
                               const bool&Protein,
                               const bool&Center,
                               const bool&Merge)const;


             void mergeInteractions();
             void loadTripletsComb(const bool full);

             static double *path_dist;
             static				 unsigned int  path_dist_size;
             static				 unsigned int  size_triplets;
             static char Int_TO_Char[10];
             static int vect_list[1000][1000];
             static bool Load_triplet;
             Fingerprint& generateTriplets(InterResults &interResult, const bool full);


             /**
              * @brief Convert a binding mode into a residu-based fingerprint
              * @param interResult : interaction set to be converted
              * @param fgpType : 0-standard IFP \
              * 1-polar IFP\
              * 2-extended IFP
              * 3-extended polar IFP
              *
              */
             void genIFP(InterResults& interResult,
                         const unsigned int& fgpType=0) const;


             bool comparaisoninterface (InterResults* i,InterResults* j) {
                 return (i->Ints.getMaxNumAtom() >j->Ints.getMaxNumAtom());
             }



/**
 * @brief Maximum distance allowed between the acceptor and the donor in an Hydrogen Bond.
 * @param dist : threshold in Angstroem
 * Default : 4.5 Angstroem
 */
inline static void  setDist_H          (const double &dist) { Dist_H      = dist;}


/**
 * @brief Maximum distance allowed between two hydrophobic atoms.
 * @param dist : in Angstroems
 * Default : 4.5 Angstroem
 */
inline static  void  setDist_Hyd       (const double &dist) { Dist_Hyd    = dist;} /*!<  */


/**
 * @brief Maximum distance allowed between an anion and a cation.
 * @param dist : in Angstroems
 * Default : 4 Angstroems
 */
inline static  void  setDist_Ionic     (const double &dist) { Dist_Ionic  = dist;} /*!<  */


/**
 * @brief Maximum distance allowed between a metal and an acceptor.
 * @param dist : in Angstroems
 *  Default : 2.8 Angstroem
 */
inline static  void  setDist_Metal     (const double &dist) { Dist_Metal  = dist;}


/**
 * @brief Maximum distance allowed between centers of two aromatic rings.
 * @param dist : in Angstroems
 *
 *  Default : 4 Angstroem
 */
inline static  void  setDist_Arom      (const double &dist) { Dist_Arom   = dist;}


/**
 * @brief Maximum distance allowed between the acceptor and the donor in an Hydrogen Bond.
 * @param dist : threshold in Angstroem
 * Default : 4.5 Angstroem
 */
inline static void  setdist_H          (const double &dist) { dist_H      = dist;}


/**
 * @brief Maximum distance allowed between two hydrophobic atoms.
 * @param dist : in Angstroems
 * Default : 4.5 Angstroem
 */
inline static  void  setdist_Hyd       (const double &dist) { dist_Hyd    = dist;} /*!<  */


/**
 * @brief Maximum distance allowed between an anion and a cation.
 * @param dist : in Angstroems
 * Default : 4 Angstroems
 */
inline static  void  setdist_Ionic     (const double &dist) { dist_Ionic  = dist;} /*!<  */


/**
 * @brief Maximum distance allowed between a metal and an acceptor.
 * @param dist : in Angstroems
 *  Default : 2.8 Angstroem
 */
inline static  void  setdist_Metal     (const double &dist) { dist_Metal  = dist;}


/**
 * @brief Maximum distance allowed between centers of two aromatic rings.
 * @param dist : in Angstroems
 *
 *  Default : 4 Angstroem
 */
inline static  void  setdist_Arom      (const double &dist) { dist_Arom   = dist;}


/**
 * @brief Maximum distance allowed between the acceptor and the donor atom in weak H-Bond
 * @param dist : in Angstroems
 *
 * Default:2.8 Angstreoms
 */
inline static void setdist_WHBond(const double& dist){Dist_WHBond=dist;}

/**
 * @brief Maximum distance allowed between the cation and the aromatic center
 * @param dist : in Angstroems
 *
 * Default: 4 Angstreoms
 */
inline static void setDist_PICation(const double& dist){Dist_PiCation = dist;}
//inline static void setdist_PICation(const double& dist){dist_PiCation=dist;}

/**
 * @brief Angle between the donor, the hydrogen and the acceptor atoms for an hydrogen bond.
 * @param angl : in radian
 *
 * Default : PI (180 degrees)
 */
inline static  void  setAngl_H         (const double &angl) { Angl_H      = angl;}


/**
 * @brief Tolerance for the hydrogen bond angle.
 * @param angl : in radian
 *
 * Default : PI/3 (60 degrees) - The range by default is 180° +/- 60
 */
inline static  void  setAngl_Tol_H     (const double &angl) { AngT_H      = angl;}


/**
 * @brief Angle between the two aromatic rings in a face to face interaction.
 * @param angl : in radian
 *
 * Default : PI (180 degrees)
 */
inline static  void  setAngl_AromFF    (const double &angl) { Angl_AromFF = angl;}


/**
 * @brief Tolerance for the face to face interaction angle. Default : PI/6 (30 degrees) - The range by default is 180 +/- 30
 * @param angl : in radian
 */
inline static void  setAngl_Tol_AromFF(const double &angl) { AngT_AromFF = angl;}


/**
 * @brief Angle between the two aromatic rings in a edge to face interaction.
 * @param angl : in radian
 *
 * Default : PI/2 (90 degrees)
 */
inline static void  setAngl_AromEF    (const double &angl) { Angl_AromEF = angl;}


/**
 * @brief Tolerance for the edge to face interaction angle.
 * @param angl : in radian
 *
 * Default : PI/3 (60 degrees) - The range by default is 90 +/- 60
 */
inline static void  setAngl_Tol_AromEF(const double &angl) { AngT_AromFF = angl;}


/**
 * @brief Angle between the vector aromatic center-cation and the normal vector of the aromatic cycle.
 * @param angl : in radian
 *
 * Default : PI (180 degrees)
 */
inline static void  setAngl_PICation(const double &angl) { Angl_PiCation = angl;}


/**
 * @brief Tolerance for the PI Cation interaction angle.
 * @param angl : in radian
 *
 * Default : PI/6 (20 degrees) - The range by default is 180 +/- 20
 */
inline static void  setAngl_Tol_PICation(const double &angl) { AngT_PiCation = angl;}



/**
 * @brief Angle between the HBond acceptor, the hydrogen and the HBond donor in Weak HBond interaction.
 * @param angl : in radian
 *
 * Default : PI (180 degrees)
 */
inline static void  setAngl_WHBond(const double &angl) { Angl_WHBond= angl;}


/**
 * @brief Tolerance for Weak HBond interaction angle
 * @param angl : in radian
 *
 * Default : PI/6 (20 degrees) - The range by default is 180 +/- 20
 */
inline static void  setAngl_Tol_WHBond(const double &angl) { AngT_WHBond = angl;}





/**
 * @brief Maximum distance allowed between the acceptor and the donor in an Hydrogen Bond.
 * @return distance in angstroems
 */
inline static  double  getDist_H         ()  { return  Dist_H      ;}


/**
 * @brief Maximum distance allowed between two hydrophobic atoms.
 * @return
 */
inline static  double  getDist_Hyd       ()  { return  Dist_Hyd    ;}


/**
 * @brief Maximum distance allowed between an anion and a cation.
 * @return
 */
inline static  double  getDist_Ionic     ()  { return  Dist_Ionic  ;}


/**
 * @brief Maximum distance allowed between a metal and an acceptor
 * @return
 */
inline static  double  getDist_Metal     ()  { return  Dist_Metal  ;}


/**
 * @brief  Maximum distance allowed between centers of two aromatic rings.
 * @return
 */
inline static  double  getDist_Arom      ()  { return  Dist_Arom   ;}


/**
 * @brief Angle between the donor, the hydrogen and the acceptor atoms for an hydrogen bond.
 * @return
 */
inline static  double  getAngl_H         ()  { return  Angl_H      ;}


/**
 * @brief Tolerance for the hydrogen bond angle.
 * @return
 */
inline static  double  getAngl_Tol_H     ()  { return  AngT_H      ;}


/**
 * @brief Angle between the two aromatic rings in a face to face interaction.
 * @return
 */
inline static  double  getAngl_AromFF    ()  { return  Angl_AromFF ;}


/**
 * @brief Tolerance for the face to face interaction angle.
 * @return
 */
inline static  double  getAngl_Tol_AromFF()  { return  AngT_AromFF ;}


/**
 * @brief Angle between the two aromatic rings in a edge to face interaction.
 * @return
 */
inline static  double  getAngl_AromEF    ()  { return  Angl_AromEF ;}


/**
 * @brief Tolerance for the edge to face interaction angle.
 * @return
 */
inline static  double  getAngl_Tol_AromEF()  { return  AngT_AromFF ;}


inline   ICMole::Grid& getGrid() { return grid;}


inline const ICMole::Coords& getCenter(const InterResults& interresult, unsigned const int& n) const {return interresult.listInters.at(n).center;}












};

}
#endif // INTERACTION_H
