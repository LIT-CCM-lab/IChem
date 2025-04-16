#ifndef MOLEDATA_H
#define MOLEDATA_H
#include <string>
#define NB_ATMS 119
#define NB_MOL2 51
#define NBAA 29
#define NB_COLOR 17
#define NB_RESTYPE 13
#define NB_MOLETYPE 18
#define NB_INTTYPE 13
#define NBHN 105
#define NB_PHARMTYPE 7

namespace MoleType {
const unsigned int
    PROTEIN=0,
    COFACTOR=1,
    LIGAND=2,
    WATER=3,
    SITE=4,
    OTHER=5,
    CAV_ALL=6,
    CAV_4=7,
    CAV_6=8,
    CAV_8=9,
    CAV_12=10,
    INT_LIG=11,
    INT_CENT=12,
    INT_PROT=13,
    INT_MERG=14,
    UNDEFINED=15,
    NUCLEIC=16,
    PHARMACOPHORE=17;
}

namespace PharType{
const unsigned int
    HBDONOR=3,
    HBACCEPTOR=2,
    IONPOS=5,
    IONNEG=7,
    AROMATIC=13,
    HYDROPHOBIC=11,
    METAL=17,
    EXCLUSION=19;
}


namespace BondType      {
/*!< This namespace constains the following constants to represents bond type : */

const unsigned int
SINGLE=101,  /*!< Single bond */
DOUBLE=102,  /*!< Double bond */
TRIPLE=103,  /*!< Triple bond */
AROMATIC=104,  /*!< Aromatic bond */
AMIDE=105,  /*!< Amide bond */
ANY=106,  /*!< Any of the above */
DUMMY=107,  /*!< Dummy bond*/
UNDEFINED=108,  /*!< Unknown */
LINKER=109,
SIDE_CHAIN=110;

}

namespace ResType {
const unsigned int

STD_AA=0,
MOD_AA=1,
NUCLEIC=2,
COFACTOR=3,
METAL=4,
ION=5,
WATER=6,
LIGAND=7,
PROSTHETIC=8,
ORGANOMET=9,
SUGAR=10,
UNWANTED=11,
UNKNOWN=12;
}

namespace FileFormat{
const unsigned int
MOL2=10001,
PDB=10002,
SDF=10003,
UNDEFINED=10004,
CHM=10005,
PML=10006
        ;
}




namespace InterType   {
const unsigned int
UNDEFINED=0,
HBOND_PROT=1,  /*!< HBond Protein donor        */
HBOND_LIG=2,  /*!< HBond Ligand donor         */
IONIC_PROT=3,  /*!< Ionic Protein cation       */
IONIC_LIG=4,  /*!< Ionic Ligand cation        */
HYDROPHOBIC=5,  /*!< Hydrophobic                */
METAL=6,  /*!< Metal acceptor interaction */
ARFACEFACE=7,  /*!< Aromatic face to face      */
AREDGEFACE=8,  /*!< Aromatic edge to face      */
PICATION=9,  /*!< PI cation                  */
WHBOND_LIG=10,
WHBOND_PROT=11,
EXCLUSION=12,   /*!< Exclusion sphere           */
METAL_ACC=13  /*!< Metal charged acceptor */
        ;
}

static const struct
{
    std::string name;
    double RED;
    double GREEN;
    double BLUE;
} ColorData[NB_COLOR]={
    {"GREEN"       ,  0,255,  0},
    {"CYAN"        ,  0,255,255},
    {"BLUE"        ,  0,  0,255},
    {"MAGENTA"     ,255,  0,255},
    {"PURPLE"      ,127,  0,255},
    {"RED"         ,255,  0,  0},
    {"ORANGE"      ,255,127,  0},
    {"YELLOW"      ,255,255,  0},
    {"WHITE"       ,255,255,255},
    {"SPRINGGREEN" ,  0,255,127},
    {"ORANGERED"   ,255, 38,  0},
    {"VIOLET"      ,170,  0,255},
    {"TAN"         ,210,180,140},
    {"GOLD"        ,255,215,  0},
    {"BLACK"       ,255,255,255},
    {"BROWN"       ,165, 42, 42},
    {"MAUVE"       ,205,130,207}

};

static const struct
{
    unsigned int num;    /*!< Atomic number */
    char nom[4];         /*!< short standard name of the atom */
    std::string longname;/*!< Standard complete name of the atom */
    unsigned int rvb[3]; /*!< Default RGB color */
    double VDWradius;    /*!< ATOMIC VdW Radii : M. Mantina, D.G. Truhlar, J. Phys. Chem. A., 2009, 113, 5806-5812 doi:10.1021/jp8111556 */
    double radius;       /*!< Dalton Trans. 10.1039/b801115j , 2008, 2832–2838 */
    double weight;       /*!< ATOMIC Weights : M.E. Wieser, Pure Appl. Chem., 2006, 78, 2051-2066 doi:10.1351/pac200678112051*/
    bool biol;          /*!< Is this atom relevant in biology*/

}

AtomData[NB_ATMS] =	{
    //  N  Symbol Name       Color:R G   B      VDW     Radius   Weight    Biol
    {1,  "H" ,"Hydrogen",    {192,192,192}, 1.50,    0.31,   1.0079,   true },
    {2,  "He","Helium",      {192,192,192}, 1.40,    0.28,   4.0026,   false}, // 0.70 taken from RDKit
    {3,  "Li","Lithium",     {192,192,192}, 1.81,    1.28,   6.9410,   false}, // Radius : 1.23 from RDKit
    {4,  "Be","Beryllium",   {192,192,192}, 1.53,    0.96,   9.0128,   false},
    {5,  "B" ,"Boron",       {192,192,192}, 1.92,    0.84,  10.8110,   false},
    {6,  "C" ,"Carbon",      { 32,128, 32}, 1.70,    0.76,  12.0107,   true },
    {7,  "N" ,"Nitrogen",    { 64, 64,255}, 1.55,    0.71,  14.0067,   true },
    {8,  "O" ,"Oxygen",      {255, 32, 32}, 1.52,    0.66,  15.9994,   true },
    {9,  "F" ,"Fluor",       {100,100,  0}, 1.47,    0.57,  18.9984,   true },
    {10, "Ne","Neon",        {100,100,  0}, 1.54,    0.58,  20.1797,   false},// 0.70 Taken from RDKit
    {11, "Na","Sodium",      {100,100,  0}, 2.27,    1.66,  32.0655,   true },
    {12, "Mg","Magnesium",   {100,100,  0}, 1.73,    1.41,  35.4532,   true }, // MODIFIED FOR VolSite : Initial 1.73
    {13, "Al","Aluminium",   {100,100,  0}, 1.84,    1.21,  26.8715,   false},
    {14, "Si","Silicon",     {100,100,  0}, 2.10,    1.11,  28.0855,   false},
    {15, "P" ,"Phosphorus",  {100,100,  0}, 1.80,    1.07,  30.9738,   true },
    {16, "S" ,"Sulfur",      {205,207,  0}, 1.80,    1.05,  32.0655,   true },
    {17, "Cl","Chloride",    {102,201,  0}, 1.75,    1.02,  35.4532,   true },
    {18, "Ar","Argon",       {102,201,  0}, 1.88,    1.06,  39.9480,   false},// 1.74 taken from RDKit
    {19, "K" ,"Potassium",   {102,201,  0}, 2.75,    2.03,  39.0983,   true },
    {20, "Ca","Calcium",     {102,201,  0}, 2.31,    1.76,  40.0784,   true },
    {21, "Sc","Scandium",    {102,201,  0}, 1.70,    1.70,  44.9559,   false},// 1.70 taken from RDKit
    {22, "Ti","Titanium",    {102,201,  0}, 1.70,    1.60,  47.8671,   false},// 1.70 taken from RDKit
    {23, "V" ,"Vanadium",    {102,201,  0}, 1.70,    1.53,  50.8415,   false},// 1.70 taken from RDKit
    {24, "Cr","Chromium",    {102,201,  0}, 1.70,    1.39,  51.9962,   false},// 1.70 taken from RDKit
    {25, "Mn","Manganese",   {102,201,  0}, 1.70,    1.39,  54.9381,   true },// 1.70 taken from RDKit
    {26, "Fe","Iron",        {102,201,  0}, 1.70,    1.32,  55.8452,   true },// 1.70 taken from RDKit
    {27, "Co","Cobalt",      {102,201,  0}, 1.70,    1.26,  58.9332,   true },// 1.70 taken from RDKit
    {28, "Ni","Nickel",      {102,201,  0}, 1.63,    1.24,  58.6934,   true },// 1.63 taken from RDKit
    {29, "Cu","Copper",      {102,201,  0}, 1.40,    1.32,  63.5463,   true },// 1.40 taken from RDKit
    {30, "Zn","Zinc",        {102,201,  0}, 1.39,    1.22,  65.4094,   true },// 1.39 taken from RDKit
    {31, "Ga","Gallium",     {102,201,  0}, 1.87,    1.22,  69.7230,   false},
    {32, "Ge","Germanium",   {102,201,  0}, 2.11,    1.20,  72.6410,   false},
    {33, "As","Arsenic",     {102,201,  0}, 1.85,    1.19,  74.9216,   false},
    {34, "Se","Selenium",    {102,201,  0}, 1.90,    1.20,  78.9600,   true},
    {35, "Br","Bromine",     {102,201,  0}, 1.83,    1.20,  79.9040,   true },
    {36, "Kr","Krypton",     {102,201,  0}, 2.02,    1.16,  83.7980,   false},
    {37, "Rb","Rubidium",    {102,201,  0}, 3.03,    2.20,  85.4678,   false},
    {38, "Sr","Strontium",   {102,201,  0}, 2.49,    1.95,  87.6210,   false},
    {39, "Y" ,"Yttrium",     {102,201,  0}, 0,       1.90,  88.9059,   false},
    {40, "Zr","Zirconium",   {102,201,  0}, 0,       1.75,  91.2242,   false},
    {41, "Nb","Niobium",     {102,201,  0}, 0,       1.64,  92.9064,   false},
    {42, "Mo","Molybdenum",  {102,201,  0}, 0,       1.54,  95.9420,   false},
    {43, "Tc","Technetium",  {102,201,  0}, 0,       1.47,   0.0000,   false},
    {44, "Ru","Ruthenium",   {102,201,  0}, 0,       1.46,  101.072,   false},
    {45, "Rh","Rhodium",     {102,201,  0}, 0,       1.42,  102.906,   false},
    {46, "Pd","Palladium",   {102,201,  0}, 0,       1.39,  106.421,   false},
    {47, "Ag","Silver",      {102,201,  0}, 0,       1.45,  107.868,   false},
    {48, "Cd","Cadmium",     {102,201,  0}, 0,       1.44,  112.412,   false},
    {49, "In","Indium",      {102,201,  0}, 1.93,    1.42,  114.818,   false},
    {50, "Sn","Tin",         {102,201,  0}, 2.17,    1.39,  118.711,   false},
    {51, "Sb","Antimony",    {102,201,  0}, 2.06,    1.39,  121.760,   false},
    {52, "Te","Tellurium",   {102,201,  0}, 2.06,    1.38,  127.603,   false},
    {53, "I" ,"Iodine",      {102,201,  0}, 1.98,    1.39,  126.904,   true },
    {54, "Xe","Xenon",       {102,201,  0}, 2.16,    1.40,  131.293,   false},
    {55, "Cs","Cesium",      {102,201,  0}, 3.43,    2.44,  132.905,   false},
    {56, "Ba","Barium",      {102,201,  0}, 2.68,    2.15,  137.327,   false},
    {57, "La","Lanthanum",   {102,201,  0}, 0,       2.07,  138.905,   false},
    {58, "Ce","Cerium",      {102,201,  0}, 0,       2.04,  140.116,   false},
    {59, "Pr","Praseodymium",{102,201,  0}, 0,       2.03,  140.907,   false},
    {60, "Nd","Neodymium",   {102,201,  0}, 0,       2.01,  144.242,   false},
    {61, "Pm","Promethium",  {102,201,  0}, 0,       1.99,  0,         false},
    {62, "Sm","Samarium",    {102,201,  0}, 0,       1.98,  150.362,   false},
    {63, "Eu","Europium",    {102,201,  0}, 0,       1.98,  151.9641,  false},
    {64, "Gd","Gadolinium",  {102,201,  0}, 0,       1.96,  157.253,   false},
    {65, "Tb","Terbium",     {102,201,  0}, 0,       1.94,  158.9253,  false},
    {66, "Dy","Dysprosium",  {102,201,  0}, 0,       1.92,  162.5001,  false},
    {67, "Ho","Holmium",     {102,201,  0}, 0,       1.92,  164.9303,  false},
    {68, "Er","Erbium",      {102,201,  0}, 0,       1.89,  167.2593,  false},
    {69, "Tm","Thulium",     {102,201,  0}, 0,       1.90,  168.9342,  false},
    {70, "Yb","Ytterbium",   {102,201,  0}, 0,       1.87,  173.043,   false},
    {71, "Lu","Lutetium",    {102,201,  0}, 0,       1.87,  174.9671,  false},
    {72, "Hf","Hafnium",     {102,201,  0}, 0,       1.75,  178.492,   false},
    {73, "Ta","Tantalum",    {102,201,  0}, 0,       1.70,  180.9479,  false},
    {74, "W" ,"Tungsten",    {102,201,  0}, 0,       1.62,  183.841,   false},
    {75, "Re","Rhenium",     {102,201,  0}, 0,       1.51,  186.2071,  false},
    {76, "Os","Osmium",      {102,201,  0}, 0,       1.44,  190.233,   false},
    {77, "Ir","Iridium",     {102,201,  0}, 0,       1.41,  192.217,   false},
    {78, "Pt","Platinum",    {102,201,  0}, 0,       1.36,  195.0849,  false},
    {79, "Au","Gold",        {102,201,  0}, 0,       1.36,  196.966,   false},
    {80, "Hg","Mercury",     {102,201,  0}, 0,       1.32,  200.592,   false},
    {81, "Tl","Thallium",    {102,201,  0}, 1.96,    1.45,  204.3833,  false},
    {82, "Pb","Lead",        {102,201,  0}, 2.02,    1.46,  207.21,    false},
    {83, "Bi","Bismuth",     {102,201,  0}, 2.07,    1.48,  208.9804,  false},
    {84, "Po","Polonium",    {102,201,  0}, 1.97,    1.40,  0,         false},
    {85, "At","Astatine",    {102,201,  0}, 2.02,    1.50,     0,         false},
    {86, "Rn","Radon",       {102,201,  0}, 2.20,    1.50,     0,         false},
    {87, "Fr","Francium",    {102,201,  0}, 3.48,    2.60,     0,         false},
    {88, "Ra","Radium"	,    {102,201,  0}, 0,2.21,226.0254,false},
    {89, "Ac","Actinium"	,{102,201,  0}, 0,2.15,227.0278,false},
    {90, "Th","Thorium"	    ,{102,201,  0}, 0,2.06,232.0381,false},
    {91, "Pa","Protactinium",{102,201,  0}, 0,2.00,231.0359,false},
    {92, "U" ,"Uranium"	    ,{102,201,  0}, 0,1.96,238.0289,false},
    {93, "Np","Neptunium"	,{102,201,  0}, 0,1.90,237.0482,false},
    {94, "Pu","Plutonium"	,{102,201,  0}, 0,1.87,244.0642,false},
    {95, "Am","Amricium"	,{102,201,  0}, 0,1.80,243.0614,false},
    {96, "Cm","Curium"	    ,{102,201,  0}, 0,1.69,247.0703,false},
    {97, "Bk","Berklium"	,{102,201,  0}, 0,0,247.0703,false},
    {98, "Cf","Californium"	,{102,201,  0}, 0,0,251.0796,false},
    {99, "Es","Einsteinium"	,{102,201,  0}, 0,0,252.0829,false},
    {100,"Fm","Fermium"	    ,{102,201,  0}, 0,0,257.0951,false},
    {101,"Md","Mendlvium"	,{102,201,  0}, 0,0,258.0986,false},
    {102,"No","Noblium"	    ,{102,201,  0}, 0,0,259.1009,false},
    {103,"Lr","Lawrencium"	,{102,201,  0}, 0,0,260.1053,false},
    {104,"Rf","Rutherfordium",{102,201,  0}, 0,0,261.1087,false},
    {105,"Db","Dubnium"	    ,{102,201,  0}, 0,0,262.1138,false},
    {106,"Sg","Seaborgium"	,{102,201,  0}, 0,0,263.1182,false},
    {107,"Bh","Bohrium"	    ,{102,201,  0}, 0,0,262.1229,false},
    {108,"Hs","Hassium"	    ,{102,201,  0}, 0,0,265 	,false},
    {109,"Mt","Meitnerium"	,{102,201,  0}, 0,0,266 	,false},
    {110,"Ds","Darmstadtium",{102,201,  0}, 0,0,269 	,false},
    {111,"Rg","Roentgenium"	,{102,201,  0}, 0,0,272 	,false},
    {112,"Cn","Copernicium"	,{102,201,  0}, 0,0,277 	,false},
    {113,"Uut","Ununtrium"	,{102,201,  0}, 0,0,0        ,false},
    {114,"Fl" ,"Flrovium"	,{102,201,  0}, 0,0,0        ,false},
    {115,"Uup","Ununpentium",{102,201,  0}, 0,0,0        ,false},
    {116,"Lv" ,"Livermorium",{102,201,  0}, 0,0,0        ,false},
    {117,"Uus","Ununseptium",{102,201,  0}, 0,0,0        ,false},
    {118,"Uuo","Ununoctium"	,{102,201,  0}, 0,0,0        ,false},
    {119,"Du" ,"Dummy"      ,{ 10, 10, 10}, 0.00,    0.00,        0.00,      true }
};


static const struct
{
    std::string mol2type;
    std::string property;    /*!<  Pharmacophoric property 1 */
    unsigned int atomic_num;  /*!< Atomic name of the mol2 type */
}
MOL2_SPEC[NB_MOL2]= {
    {"H",    "NOP",          1},
    {"C.3",  "H, np, ",      6},    //
    {"C.2",  "np",           6},    //
    {"O.2",  "A",            8},//
    {"N.am", "D",            7},//
    {"C.ar", "Ar, np, ",     6},
    {"O.3",  "A, D, ",       8},//
    {"O.co2","A, MOINS, ",   8},//
    {"N.pl3","D",            7},//
    {"S.3",  "H, A, np, ",   16},    //
    {"N.4",  "D, PLUS, ",    7},//
    {"N.2",  "A",            7},//
    {"C.cat","PLUS, np, ",         6},//
    {"P.3",  "NOP",         15},
    {"N.ar", "Ar, A, D, ",   7},
    {"F",    "H",            9},//
    {"Mg",   "D, PLUS, ",   12},//
    {"S.o2", "NOP",      16},
    {"Cl",   "H, np, ",     17},//
    {"I",    "H, np, ",      53},//
    {"Zn",   "D, PLUS, ",   30},//
    {"N.3",  "D",            7},//
    {"C.1",  "H, np, ",      6},//
    {"Mn",   "D, PLUS, ",   25},//
    {"K",    "",            19},
    {"Br",   "H, np, ",    35},//
    {"Ca",   "D, PLUS, ",   20},//
    {"S.2",  "",           16},
    {"N.1",  "A",            7},//
    {"Fe",   "D, PLUS, ",   26},//
    {"S.o",  "NOP",         16},
    {"Du",   "",             1},
    {"S",    "" ,           16},//
    {"P",    "",            15},
    {"Na",   "PLUS",        11},
    {"Li",   "",             3},
    {"Al",   "",            13},
    {"Cu",   "D, PLUS, ",   29},//
    {"LP",   "",             1},
    {"Du.C", "",             1},
    {"O.spc","A, D, ",       8},//
    {"O.t3p","A, D, ",       8},//
    {"H.spc","NOP",          1},
    {"H.t3p","NOP",          1},
    {"Si",   "",            14},
    {"Cr.oh","",            24},
    {"Cr.th","",            24},
    {"Se",   "",            34},
    {"Sn",   "",            50},
    {"Mo",   "",            42},
    {"Co.oh","",            27},

};





static const struct
{
    std::string code;       /*!< Amino Acid Code */
    std::string name;/* Full name of the amino acid*/
    double MW;              /* Molecular weight of AA */
    std::string property;
    unsigned int nb_heavy_atom;
    unsigned int nb_atoms;
    unsigned int restype;

}
AAcid[NBAA]={
    {"A", "ALA", 89.0940, "AP, AL", 5,10, ResType::STD_AA},
    {"C", "CYS",121.1540, "PO"    , 7,14, ResType::STD_AA},
    {"D", "ASP",133.1038, "AC"    , 9,16, ResType::STD_AA},
    {"E", "GLU",147.1307, "AC"    ,10,19, ResType::STD_AA},
    {"F", "PHE",165.1918, "AP, AR",12,23, ResType::STD_AA},
    {"G", "GLY", 75.0671, "AP"    , 5,10, ResType::STD_AA},
    {"H", "HIS",155.1563, "BA, AR",11,20, ResType::STD_AA},
    {"I", "ILE",131.1746, "AP, AL", 9,22, ResType::STD_AA},
    {"K", "LYS",146.1893, "BA"    ,10,24, ResType::STD_AA},
    {"L", "LEU",131.1746, "AP, AL", 9,22, ResType::STD_AA},
    {"M", "MET",149.2078, "AP"    , 9,20, ResType::STD_AA},
    {"N", "ASN",132.1190, "PO"    , 9,17, ResType::STD_AA},
    {"O", "PYL",255.3134, "PO"    ,18,39, ResType::STD_AA},
    {"P", "PRO",115.1305, ""      , 8,17, ResType::STD_AA},
    {"Q", "GLN",146.1459, "PO"    ,10,20, ResType::STD_AA},
    {"R", "ARG",174.2027, "BA"    ,12,26, ResType::STD_AA},
    {"S", "SER",105.0934, "PO"    , 7,14, ResType::STD_AA},
    {"T", "THR",119.1203, "PO"    , 8,17, ResType::STD_AA},
    {"U", "SEC",168.0530, "PO"    , 8,15, ResType::STD_AA},
    {"V", "VAL",117.1478, "AP, AL", 8,19, ResType::STD_AA},
    {"W", "TRP",204.2284, "AP, AR",14,26, ResType::STD_AA},
    {"Y", "TYR",181.1910, "PO, AR",13,24, ResType::STD_AA},
    {"" , "HIE",       0, ""      ,0 ,0 , ResType::STD_AA},
    {"" , "HID",       0, ""      ,0 ,0 , ResType::STD_AA},
    {"" , "HOH",   18.01, ""      ,1, 3 , ResType::WATER },
    {"" , "WAT",   18.01, ""      ,1, 3 , ResType::WATER },
    {"" , "WTR",   18.01, ""      ,1, 3,  ResType::WATER },
    {"" , "MTO",   18.01, ""      ,1, 3,  ResType::WATER },
    {"" , "DOD",   18.01, ""      ,1, 3,  ResType::WATER },
};


#ifdef GUI_DESIGN
static const struct
{
    std::string name;
    double R;
    double G;
    double B;
} ColorData[NB_COLOR]={
    {"GREEN"       ,  0,255,  0},
    {"CYAN"        ,  0,255,255},
    {"BLUE"        ,  0,  0,255},
    {"MAGENTA"     ,255,  0,255},
    {"PURPLE"      ,127,  0,255},
    {"RED"         ,255,  0,  0},
    {"ORANGE"      ,255,127,  0},
    {"YELLOW"      ,255,255,  0},
    {"WHITE"       ,255,255,255},
    {"SPRINGGREEN" ,  0,255,127},
    {"ORANGERED"   ,255, 38,  0},
    {"VIOLET"      ,170,  0,255},
    {"TAN"         ,210,180,140},
    {"GOLD"        ,255,215,  0},
    {"BLACK"       ,255,255,255},
    {"BROWN"       ,165, 42, 42},
    {"MAUVE"       ,205,130,207}

};
#endif



static const struct
{
    std::string atmName;       /*!< Amino Acid Code */
    std::string resName    ;/* Full name of the amino acid*/
    std::string H1;
    std::string H2;
    std::string H3;
}
AtomName[NBHN]=
{
    {"CA" ,"ALA","HA"  ,""    ,""    },
    {"CB" ,"ALA","HB1" ,"HB2" ,"HB3" },
    {"CA" ,"ARG","HA"  ,""    ,""    },
    {"N"  ,"ARG","HNCA",""    ,""    },
    {"CB" ,"ARG","HB2" ,"HB3" ,""    },
    {"CG" ,"ARG","HG2" ,"HG3" ,""    },
    {"CD" ,"ARG","HD2" ,"HD3" ,""    },
    {"NE" ,"ARG","HE"  ,""    ,""    },
    {"NH1","ARG","HH11","HH12",""    },
    {"NH2","ARG","HH21","HH22",""    },
    {"CA" ,"ASN","HA"  ,""    ,""    },
    {"N"  ,"ASN","HNCA",""    ,""    },
    {"CB" ,"ASN","HB2" ,"HB3" ,""    },
    {"ND2","ASN","HD21","HD22",""    },
    {"CA" ,"ASP","HA"  ,""    ,""    },
    {"N"  ,"ASP","HNCA",""    ,""    },
    {"CB" ,"ASP","HB2" ,"HB3" ,""    },
    {"CA" ,"CYS","HA"  ,""    ,""    },
    {"N"  ,"CYS","HNCA",""    ,""    },
    {"CB" ,"CYS","HB2" ,"HB3" ,""    },
    {"SG" ,"CYS","HG"  ,""    ,""    },
    {"CA" ,"GLU","HA"  ,""    ,""    },
    {"N"  ,"GLU","HNCA",""    ,""    },
    {"CB" ,"GLU","HB2" ,"HB3" ,""    },
    {"CG" ,"GLU","HG2" ,"HG3" ,""    },
    {"CA" ,"GLN","HA"  ,""    ,""    },
    {"N"  ,"GLN","HNCA",""    ,""    },
    {"CB" ,"GLN","HB2" ,"HB3" ,""    },
    {"CG" ,"GLN","HG2" ,"HG3" ,""    },
    {"NE2","GLN","HE21","HE22",""    },
    {"CA" ,"GLY","HA2" ,"HA3" ,""    },
    {"N"  ,"GLY","HNCA",""    ,""    },
    {"CA" ,"HIS","HA"  ,""    ,""    },
    {"N"  ,"HIS","HNCA",""    ,""    },
    {"CB" ,"HIS","HB2" ,"HB3" ,""    },
    {"ND1","HIS","HD1" ,""    ,""    },
    {"CE1","HIS","HE1" ,""    ,""    },
    {"CD2","HIS","HD2" ,""    ,""    },
    {"CA" ,"ILE","HA"  ,""    ,""    },
    {"N"  ,"ILE","HNCA",""    ,""    },
    {"CB" ,"ILE","HB"  ,""    ,""    },
    {"CG2","ILE","HG21","HG22","HG23"},
    {"CG1","ILE","HG12","HG13",""    },
    {"CD1","ILE","HD11","HD12","HD13"},
    {"CA" ,"LEU","HA"  ,""    ,""    },
    {"N"  ,"LEU","HNCA",""    ,""    },
    {"CB" ,"LEU","HB2" ,"HB3" ,""    },
    {"CG" ,"LEU","HG"  ,""    ,""    },
    {"CD2","LEU","HD21","HD22","HD23"},
    {"CD1","LEU","HD11","HD12","HD13"},
    {"CA" ,"LYS","HA"  ,""    ,""    },
    {"N"  ,"LYS","HNCA",""    ,""    },
    {"CB" ,"LYS","HB2" ,"HB3" ,""    },
    {"CG" ,"LYS","HG2" ,"HG3" ,""    },
    {"CD" ,"LYS","HD2" ,"HD3" ,""    },
    {"CE" ,"LYS","HE2" ,"HE3" ,""    },
    {"NZ" ,"LYS","HZ1" ,"HZ2" ,"HZ3" },
    {"CA" ,"MET","HA"  ,""    ,""    },
    {"N"  ,"MET","HNCA",""    ,""    },
    {"CB" ,"MET","HB2" ,"HB3" ,""    },
    {"CG" ,"MET","HG2" ,"HG3" ,""    },
    {"CE" ,"MET","HE1" ,"HE2" ,"HE3" },
    {"CA" ,"PHE","HA"  ,""    ,""    },
    {"N"  ,"PHE","HNCA",""    ,""    },
    {"CB" ,"PHE","HB2" ,"HB3" ,""    },
    {"CD2","PHE","HD2" ,""    ,""    },
    {"CE2","PHE","HE2" ,""    ,""    },
    {"CZ" ,"PHE","HZ"  ,""    ,""    },
    {"CE1","PHE","HE1" ,""    ,""    },
    {"CD1","PHE","HD1" ,""    ,""    },
    {"CA" ,"PRO","HA"  ,""    ,""    },
    {"CD" ,"PRO","HD2" ,"HD3" ,""    },
    {"CG" ,"PRO","HG2" ,"HG3" ,""    },
    {"CB" ,"PRO","HB2" ,"HB3" ,""    },
    {"CA" ,"SER","HA"  ,""    ,""    },
    {"N"  ,"SER","HNCA",""    ,""    },
    {"CB" ,"SER","HB2" ,"HB3" ,""    },
    {"OG" ,"SER","HG"  ,""    ,""    },
    {"CA" ,"THR","HA"  ,""    ,""    },
    {"N"  ,"THR","HNCA",""    ,""    },
    {"CB" ,"THR","HB"  ,""    ,""    },
    {"CG2","THR","HG21","HG22","HG23"},
    {"OG1","THR","HG1" ,""    ,""    },
    {"CA" ,"TRP","HA"  ,""    ,""    },
    {"N"  ,"TRP","HNCA",""    ,""    },
    {"CB" ,"TRP","HB2" ,"HB3" ,""    },
    {"CD1","TRP","HD1" ,""    ,""    },
    {"NE1","TRP","HE1" ,""    ,""    },
    {"CZ2","TRP","HZ2" ,""    ,""    },
    {"CH2","TRP","HH2" ,""    ,""    },
    {"CZ3","TRP","HZ3" ,""    ,""    },
    {"CE3","TRP","HE3" ,""    ,""    },
    {"CA" ,"TYR","HA"  ,""    ,""    },
    {"N"  ,"TYR","HNCA",""    ,""    },
    {"CB" ,"TYR","HB2" ,"HB3" ,""    },
    {"CD1","TYR","HD1" ,""    ,""    },
    {"CE1","TYR","HE1" ,""    ,""    },
    {"CD2","TYR","HD2" ,""    ,""    },
    {"CE2","TYR","HE2" ,""    ,""    },
    {"OH" ,"TYR","HH"  ,""    ,""    },
    {"CA" ,"VAL","HA"  ,""    ,""    },
    {"N"  ,"VAL","HNCA",""    ,""    },
    {"CB" ,"VAL","HB"  ,""    ,""    },
    {"CG2","VAL","HG21","HG22","HG23"},
    {"C12","VAL","HG11","HG12","HG13"}
};
#endif // MOLEDATA_H
