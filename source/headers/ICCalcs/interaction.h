#ifndef INTERACTION_H
#define INTERACTION_H


#include "headers/ICMole/molecule.h"
#include "headers/ICMole/fingerprint.h"
#include "headers/ICCalcs/neighborsearch.hpp"


namespace ICMole {

  struct InteractionParameters {

    // Distances max
    double Dist_H = 3.5; // HBond maximal distance
    double Dist_Hyd = 4.5; // Hydrophobic maximal distance
    double Dist_Ionic = 4.0; // Ionic maximal distance
    double Dist_Metal = 2.8; // Metal/Acceptor maximal distance
    double Dist_Arom = 5.0; // Aromatic maximal distance between center
    double Dist_PiCation = 5.0; // PI-Cation maximal distance
    double Dist_WHBond = 3.5; // Weak H-Bond maximal distance
    
    // Distances min
    double dist_H = 2.3; // HBond minimal distance
    double dist_Hyd = 3.2; // Hydrophobic minimal distance
    double dist_Ionic = 2.3; // Ionic minimal distance
    double dist_Metal = 1.8; // Metal/Acceptor minimal distance
    double dist_Arom = 3.2; // Aromatic minimal distance between center
    
    
    // Angles
    double Angl_H = M_PI; // Median angle for HBond interactions
    double AngT_H = M_PI / 3; // +/- possible deviation for HBond angle
    double Angl_AromFF = M_PI; // Median Aromatic Face to Face angle
    double AngT_AromFF = M_PI / 6; // +/- possible deviation for Aromatic Face to Face angle
    double Angl_AromEF = M_PI / 2; // Median Aromatic Edge to Face angle
    double AngT_AromEF = M_PI / 3; // possible deviation for Aromatic Edge to Face angle
    double Angl_PiCation = M_PI; // Median angle for PI-Cation interaction
    double AngT_PiCation = M_PI / 6; // Possible deviation for PI-Cation interaction
    double Angl_WHBond = M_PI; // Median angle for Weak H-Bond interaction
    double AngT_WHBond = M_PI / 6; // Possible deviation for Weak H-Bond interaction
  };

  struct resbest {
    Atom *atmP = nullptr;
    Atom *atmL = nullptr;
    double dist = std::numeric_limits<double>::max();
    unsigned short id = 0;
  };


  struct InterPoint {
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


      InterPoint(const int& p, Atom* const protein, Atom* const ligand, const Coords& center, const unsigned int& int_type, const double& dist, const double& angl= -100000) :
        point(p),
        Prot_Ref(protein != nullptr ? protein : nullptr),
        Lig_Ref(ligand != nullptr ? ligand : nullptr),   
        center(center),
        interaction(int_type),
        dist(dist),
        angle(angl),                
        merged_to(-1) {}   

  };


  struct InterResults {
    Molecule Ints;
    unsigned short CA, CZ, O, N, NZ, OD1, Zn;
    std::vector<InterPoint> listInters;
    std::string IFPString;
    Fingerprint IFP;
    std::string chainA;
    std::string chainB;

    InterResults():CA(0),CZ(0),O(0),N(0),NZ(0),OD1(0),Zn(0),chainA(""),chainB(""){}
  };


  class Interactions {

    Complex& complex;
    // Grid &grid;

    InteractionParameters params;


      /**
      * @brief Filtering rules to detect only some interactions
      * Set to true when the interaction needed to be detected
      */
      bool wInterType[NB_INTTYPE];


    /**
     * @brief Merge
     * hydrophobic interactions
     * @param interResult : interaction set to be merged
     */
    void mergeInteractions(InterResults& interResult)const;


  public:
      /**
      * @brief Constructor
      * @param complex with a least one molecule set a protein
      */
      Interactions(ICMole::Complex &complex);


      /**
      * @brief Determine interaction molecule buriedness
      */
      void calcenfouiss(std::string mol_name);


      /**
      * @brief Adds a detected interaction to the interaction result container
      *
      * Creates a new interaction entry between ligand and protein atoms, storing the
      * interaction type, distance, and angle if there is. The interaction center is
      * defined as the midpoint between the two atom coordinates. Also updates the
      * corresponding interaction counters based on interaction type.
      *
      * @param interResult Interaction result structure where the new interaction is added
      * @param atomP Protein atom involved in the interaction
      * @param atomL Ligand atom involved in the interaction
      * @param dist  Distance between the two atoms
      * @param NInter Incremental interaction index counter
      * @param angle Optional pointer to the angle value; if null, angle is set to -100000]
      * @param interactionType Type of interaction (e.g. hydrogen bond, metal, ionic, etc)
      */
      void addInteraction(InterResults& interResult, Atom& atomP, Atom& atomL, double dist, int& NInter, double* angle, unsigned int interactionType) const;
      
      /**
      * @brief Processes and filters hydrophobic interactions between ligand and protein atoms
      *
      * Identifies hydrophobic contacts between nonpolar ligand and protein atoms
      * Applies distance constraints and optionally evaluates the local hydrophobic
      * environment using a KD-tree search to ensure a sufficiently hydrophobic
      * protein neighborhood
      *
      * @param atomL Ligand atom being is nonpolar
      * @param atomP Protein atom being tested nonpolar
      * @param dist  Distance between ligand and protein atoms.
      * @param oldh  If false, local hydrophobic density is evaluated around the protein atom
      * @param interResult Interaction result container where detected interactions are stored
      * @param hydlist Map of residues to their best hydrophobic contact, used to merge close contacts
      * @param neighborSearch KD-tree object used for spatial neighborhood queries
      * @param proteinAtoms List of protein atoms used in neighbor search results
      */
      void processHydrophobicInteraction(Atom& atomL, Atom& atomP, double dist, bool oldh, InterResults& interResult, std::map<Residu*, ICMole::resbest>& hydlist, NeighborSearch& neighborSearch, const std::vector<Atom*>& proteinAtoms) const;        
      
      /**
      * @brief Detects hydrophobic interactions between aromatic rings of ligand and protein
      *
      * Compares all aromatic cycles in the ligand and protein, checking for hydrophobic
      * contacts between their atoms within the allowed distance range
      * Each valid atom pair between two aromatic rings is recorded as a hydrophobic interaction
      *
      * @param ligand Ligand molecule containing potential aromatic rings
      * @param protein Protein molecule containing potential aromatic rings
      * @param interResult Container where detected interactions are stored
      * @param NInter Incremental interaction index counter
      * @param dist_H Minimum allowed hydrophobic distance
      * @param Dist_H Maximum allowed hydrophobic distance
      */
      void checkAromaticHydrophobicInteractions(Molecule& ligand, Molecule& protein, InterResults& interResult, int& NInter, double dist_H, double Dist_H) const;
    
      /**
      * @brief Detects aromatic-aromatic and pi-cation interactions between ligand and protein rings ( I already have the function, this part must be deleted)
      *
      * Uses a KD-tree spatial neighbor search around ligand aromatic ring centers
      * to identify nearby protein aromatic rings or cationic atoms. Depending on
      * geometry and orientation, records either face-to-face, edge-to-face, or pi-cation
      * interactions.
      *
      * @param ligand Ligand molecule containing aromatic rings
      * @param neighborSearch KD-tree used to query protein atoms near ligand ring centers
      * @param proteinAtoms List of protein atoms indexed by the KD-tree
      * @param max_allowed_dist Maximum allowed distance between ring centers
      * @param interResult Interaction result container where interactions are stored
      * @param wInterType Boolean array controlling which interaction types are enabled
      * @param NInter Incremental interaction index counter
      * @param min_allowed_dist Minimum allowed distance between ring centers
      */
      void processAromaticInteractions(Molecule& ligand, NeighborSearch& neighborSearch, const std::vector<Atom*>& proteinAtoms, double max_allowed_dist, InterResults& interResult, const bool* wInterType, int& NInter, double min_allowed_dist) const;

     
      /**
      * @brief Checks for a case of metal interaction involving nitrogen and link to sulfonamide atoms
      *
      * Detects metal interactions between ligand atoms (nitrogen, sulfonamide)
      * and metal ions in the protein (e.g Zn, Fe, Mg) and applies metal specific distance thresholds
      * and records valid interactions.
      *
      * @param atomL Ligand atom to evaluate a potential nitrogen donor
      * @param atomP Protein atom is a metal
      * @param dist  Distance between the two atoms
      * @param interResult Interaction container to store results
      * @param NInter Interaction index counter
      */
      void checkMetalNitrogenSulfonamideCase(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;

      /**
      * @brief Detects general metal interactions between ligand and protein atoms
      *
      * This function checks whether a ligand atom is an acceptor and the protein atom is a metal and if these two are within the distance range
      * then it records a metal interaction. Also When the distance is slightly below 3.4 and the atomic types suggest possible coordination 
      * (e.g. metal acceptor, metal anion, or nitrogen near a metallic atom), we set the property metal acceptor to true
      *
      * @param atomL Ligand atom being tested as an acceptor
      * @param atomP Protein atom acts as a metal
      * @param dist  Distance between ligand and metal atom
      * @param interResult Container where the interaction will be recorded
      * @param NInter Running interaction counter
      */
      void checkMetalInteractions(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;

      /**
      * @brief Detects hydrogen bonds where the ligand atom acts as an acceptor and the protein atom as a donor
      *
      * Applies geometric criteria (distance and angle) to identify a hydrogen bond
      *
      * @param atomL Ligand atom acting as hydrogen bond acceptor
      * @param atomP Protein atom acting as hydrogen bond donor
      * @param dist  Distance between donor and acceptor
      * @param interResult Interaction container where results are stored
      * @param NInter Incremental interaction counter
      */
      void checkHydrogenBondLigandAcceptor(Atom& atomL, Atom& atomP, double dist,  InterResults& interResult, int& NInter) const;

      /**
      * @brief Detects hydrogen bonds where the ligand atom acts as a donor and the protein atom acts as an acceptor
      *
      * Applies geometric criteria (distance and angle) to identify a hydrogen bond
      *
      * @param atomL Ligand atom being tested as hydrogen donor
      * @param atomP Protein atom being tested as acceptor
      * @param dist  Distance between the two atoms
      * @param interResult Interaction container where results are stored
      * @param NInter Incremental interaction index counter
      */
      void checkHydrogenBondLigandDonor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;



      /**
       * @brief Detects hydrogen bonds where the ligand atom acts as a acceptor and the protein atom acts as a weak donor
       *
       * Applies geometric criteria (distance and angle) to identify a weak hydrogen bond
       *
       * @param atomL Ligand atom being tested as hydrogen acceptor
       * @param atomP Protein atom being tested as weak donor
       * @param dist  Distance between the two atoms
       * @param interResult Interaction container where results are stored
       * @param NInter Incremental interaction index counter
       */
      void checkWeakHydrogenBondLigandAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;
      
      
      /**
       * @brief Detects weak hydrogen bonds where the ligand atom acts as a weak acceptor and the protein atom acts as a weak donor or a donor
       *
       * Applies geometric criteria (distance and angle) to identify a weak hydrogen bond
       *
       * @param atomL Ligand atom being tested as a weak acceptor
       * @param atomP Protein atom being tested as weak donor or a donor
       * @param dist  Distance between the two atoms
       * @param interResult Interaction container where results are stored
       * @param NInter Incremental interaction index counter
       */
      void checkWeakHydrogenBondLigandWeakAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;

      /**
       * @brief Detects weak hydrogen bonds where the ligand atom acts as a donor and the protein atom acts as a weak acceptor
       *
       * Applies geometric criteria (distance and angle) to identify a weak hydrogen bond
       *
       * @param atomL Ligand atom being tested as hydrogen donor
       * @param atomP Protein atom being tested as a weak acceptor
       * @param dist  Distance between the two atoms
       * @param interResult Interaction container where results are stored
       * @param NInter Incremental interaction index counter
       */
      void checkWeakHydrogenBondLigandDonorProteinWeakAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;
      
      /**
       * @brief Detects weak hydrogen bonds where the ligand atom acts as a weak donor and the protein atom acts as an acceptor or a weak acceptor
       *
       * Applies geometric criteria (distance and angle) to identify a weak hydrogen bond
       *
       * @param atomL Ligand atom being tested as weak hydrogen donor
       * @param atomP Protein atom being tested as weak acceptor
       * @param dist  Distance between the two atoms
       * @param interResult Interaction container where results are stored
       * @param NInter Incremental interaction index counter
       */
      void checkWeakHydrogenBondLigandWeakDonorProteinAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;

      
      /**
      * @brief Detects ionic interactions where the ligand atom is anion and the protein is either cation or metal 
      *
      * The geometric criteria distance must be in range to detect the intertaction
      *
      * @param atomL Ligand atom acts as an anion
      * @param atomP Protein atom acts as either cation or metal 
      * @param dist  Distance between the two atoms
      * @param interResult Interaction container where results are stored
      * @param NInter Incremental interaction index counter
      */
      void checkIonicProteinInteractions(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;
      
      /**
      * @brief Detects ionic interactions where the ligand atom is cation and the protein is anion
      *
      * The geometric criteria distance must be in range to detect the intertaction
      *
      * @param atomL Ligand atom acts as a cation
      * @param atomP Protein atom acts as an anion 
      * @param dist  Distance between the two atoms
      * @param interResult Interaction container where results are stored
      * @param NInter Incremental interaction index counter
      */
      void checkIonicLigandInteractions(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;


      /**
      * @brief Detects Pi-cation interactions between ligand cationic atoms and aromatic rings in the protein
      *
      * Applies geometric and distance criteria to identify pi-cation interactions,
      * where a positively charged ligand atom is positioned near the center of
      * an aromatic ring (represented by the dummy atom "DuCy" in the protein)
      * The angle between the aromatic ring and the ligand ring vector
      * must fall within the allowed tolerance range for this interaction
      *
      * @param atomL Ligand atom being tested as a cation
      * @param atomP Protein dummy atom representing the aromatic ring center
      * @param dist  Distance between the ligand atom and the ring center
      * @param interResult Interaction container where detected interactions are stored
      * @param NInter Incremental interaction index counter
      */
      void checkPiCationInteraction(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const;


      /**
      * @brief Detect interaction according to geometric rules
      * @param ligand : molecule to detect interaction with
      * @param interResult : Detected interactions are score here
      * @param wMerge : Merge hydrophobic interactions
      */
      void detectInteractions(Molecule& ligand, InterResults& interResult, bool wMerge = true, bool oldh = true) const;


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
      void interToMOL2(  InterResults& interResult, const bool&Ligand, const bool&Protein, const bool&Center, const bool&Merge)const;

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
      void genIFP(InterResults& interResult, const unsigned int& fgpType=0) const;


      bool comparaisoninterface (InterResults* i,InterResults* j) {
          return (i->Ints.getMaxNumAtom() >j->Ints.getMaxNumAtom());
      }


      /**
       * @brief Maximum distance allowed between the acceptor and the donor in an Hydrogen Bond.
       * @param dist : threshold in Angstroem
       * Default : 4.5 Angstroem
      */
      void setDist_H(double dist) { params.Dist_H = dist; } 


      /**
       * @brief Maximum distance allowed between two hydrophobic atoms.
       * @param dist : in Angstroems
       * Default : 4.5 Angstroem
       */
      void setDist_Hyd(double dist) { params.Dist_Hyd = dist; } /*!<  */


      /**
       * @brief Maximum distance allowed between an anion and a cation.
       * @param dist : in Angstroems
       * Default : 4 Angstroems
       */
      void setDist_Ionic(double dist) { params.Dist_Ionic = dist; } /*!<  */


      /**
       * @brief Maximum distance allowed between a metal and an acceptor.
       * @param dist : in Angstroems
       *  Default : 2.8 Angstroem
       */
      void setDist_Metal(double dist) { params.Dist_Metal = dist; }


      /**
       * @brief Maximum distance allowed between centers of two aromatic rings.
       * @param dist : in Angstroems
       *  Default : 4 Angstroem
       */
      void setDist_Arom(double dist) { params.Dist_Arom = dist; }


      /**
       * @brief Maximum distance allowed between the acceptor and the donor in an Hydrogen Bond.
       * @param dist : threshold in Angstroem
       * Default : 4.5 Angstroem
       */
      void setdist_H(double dist) { params.dist_H = dist; }


      /**
       * @brief Maximum distance allowed between two hydrophobic atoms.
       * @param dist : in Angstroems
       * Default : 4.5 Angstroem
       */
      void setdist_Hyd(double dist) { params.dist_Hyd = dist; } /*!<  */


      /**
       * @brief Maximum distance allowed between an anion and a cation.
       * @param dist : in Angstroems
       * Default : 4 Angstroems
       */
      void setdist_Ionic(double dist) { params.dist_Ionic = dist; } /*!<  */


      /**
       * @brief Maximum distance allowed between a metal and an acceptor.
       * @param dist : in Angstroems
       *  Default : 2.8 Angstroem
       */
      void setdist_Metal(double dist) { params.dist_Metal = dist; }


      /**
       * @brief Maximum distance allowed between centers of two aromatic rings.
       * @param dist : in Angstroems
       *
       *  Default : 4 Angstroem
       */
      void setdist_Arom(double dist) { params.dist_Arom = dist; }


      /**
       * @brief Maximum distance allowed between the acceptor and the donor atom in weak H-Bond
       * @param dist : in Angstroems
       *
       * Default: 2.8 Angstreoms
       */
      void setdist_WHBond(double dist){ params.Dist_WHBond = dist; }

      /**
       * @brief Maximum distance allowed between the cation and the aromatic center
       * @param dist : in Angstroems
       *
       * Default: 4 Angstreoms
       */
      void setDist_PICation(double dist){ params.Dist_PiCation = dist; }

      /**
       * @brief Angle between the donor, the hydrogen and the acceptor atoms for an hydrogen bond.
       * @param angl : in radian
       *
       * Default : PI (180 degrees)
       */
      void setAngl_H(double angl) { params.Angl_H = angl; }


      /**
       * @brief Tolerance for the hydrogen bond angle.
       * @param angl : in radian
       *
       * Default : PI/3 (60 degrees) - The range by default is 180° +/- 60
       */
      void setAngl_Tol_H(double angl) { params.AngT_H = angl; }


      /**
       * @brief Angle between the two aromatic rings in a face to face interaction.
       * @param angl : in radian
       *
       * Default : PI (180 degrees)
       */
      void setAngl_AromFF(double angl) { params.Angl_AromFF = angl; }


      /**
       * @brief Tolerance for the face to face interaction angle. Default : PI/6 (30 degrees) - The range by default is 180 +/- 30
       * @param angl : in radian
       */
      void setAngl_Tol_AromFF(double angl) { params.AngT_AromFF = angl; }


      /**
       * @brief Angle between the two aromatic rings in a edge to face interaction.
       * @param angl : in radian
       *
       * Default : PI/2 (90 degrees)
       */
      void setAngl_AromEF(double angl) { params.Angl_AromEF = angl; }


      /**
       * @brief Tolerance for the edge to face interaction angle.
       * @param angl : in radian
       *
       * Default : PI/3 (60 degrees) - The range by default is 90 +/- 60
       */
      void setAngl_Tol_AromEF(double angl) { params.AngT_AromFF = angl; }


      /**
       * @brief Angle between the vector aromatic center-cation and the normal vector of the aromatic cycle.
       * @param angl : in radian
       *
       * Default : PI (180 degrees)
       */
      void setAngl_PICation(double angl) { params.Angl_PiCation = angl; }


      /**
       * @brief Tolerance for the PI Cation interaction angle.
       * @param angl : in radian
       *
       * Default : PI/6 (20 degrees) - The range by default is 180 +/- 20
       */
      void setAngl_Tol_PICation(double angl) { params.AngT_PiCation = angl; }



      /**
       * @brief Angle between the HBond acceptor, the hydrogen and the HBond donor in Weak HBond interaction.
       * @param angl : in radian
       *
       * Default : PI (180 degrees)
       */
      void setAngl_WHBond(double angl) { params.Angl_WHBond = angl; }


      /**
       * @brief Tolerance for Weak HBond interaction angle
       * @param angl : in radian
       *
       * Default : PI/6 (20 degrees) - The range by default is 180 +/- 20
       */
      void  setAngl_Tol_WHBond(double angl) { params.AngT_WHBond = angl; }


      /**
       * @brief Maximum distance allowed between the acceptor and the donor in an Hydrogen Bond.
       * @return distance in angstroems
       */
        double getDist_H() const { return  params.Dist_H; }


      /**
       * @brief Maximum distance allowed between two hydrophobic atoms.
       * @return
       */
      double getDist_Hyd() const { return  params.Dist_Hyd; }


      /**
       * @brief Maximum distance allowed between an anion and a cation.
       * @return
       */
      double getDist_Ionic() const  { return params.Dist_Ionic; }


      /**
       * @brief Maximum distance allowed between a metal and an acceptor
       * @return
       */
      double getDist_Metal() const  { return params.Dist_Metal; }


      /**
       * @brief  Maximum distance allowed between centers of two aromatic rings.
       * @return
       */
      double getDist_Arom() const  { return params.Dist_Arom; }


      /**
       * @brief Angle between the donor, the hydrogen and the acceptor atoms for an hydrogen bond.
       * @return
       */
      double getAngl_H() const { return params.Angl_H; }


      /**
       * @brief Tolerance for the hydrogen bond angle.
       * @return
       */
      double getAngl_Tol_H() const { return params.AngT_H; }


      /**
       * @brief Angle between the two aromatic rings in a face to face interaction.
       * @return
       */
      double getAngl_AromFF() const { return params.Angl_AromFF ;}


      /**
       * @brief Tolerance for the face to face interaction angle.
       * @return
       */
      double getAngl_Tol_AromFF() const { return params.AngT_AromFF; }


      /**
       * @brief Angle between the two aromatic rings in a edge to face interaction.
       * @return
       */
      double getAngl_AromEF() const { return params.Angl_AromEF; }


      /**
       * @brief Tolerance for the edge to face interaction angle.
       * @return
       */
      double getAngl_Tol_AromEF() const { return params.AngT_AromFF; }


      // inline ICMole::Grid& getGrid() { return grid; }
      inline const ICMole::Coords& getCenter(const InterResults& interresult, unsigned const int& n) const {return interresult.listInters.at(n).center;}

  };
}

#endif // INTERACTION_H
