#include "headers/ICMole/fingerprint.h"


using namespace std;
using namespace ICMole;
unsigned int Fingerprint::debug=0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// CONSTRUCTORS ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Fingerprint::Fingerprint():numeric(false)
{

}
Fingerprint::Fingerprint(const unsigned int& size,const bool& numerical):bit_size(size),name("")
{
    numeric=numerical;
    bitstring.reserve(size);
    for (unsigned int Iv = 0; Iv < size; Iv++) bitstring.push_back(0);
}

Fingerprint::Fingerprint(const string& Fname,const bool& numerical) throw(MoleExcept)
{
    bit_size=0;
    string ligne;
    numeric=numerical;
    ifstream finput;
    finput.open(Fname.c_str(),ios::in);
    if (finput.is_open() == false)  { throw MoleExcept(1111001,"Fingerprint::Constructor","Cannot open fingerprint  file : " + Fname);}
    std::getline( finput, ligne );
    size_t pos=ligne.find_first_of(":");
    if (pos != string::npos) {finput.close();loadSVM_File(Fname);return;}
    pos = ligne.find_first_of("[");
    if (pos != string::npos) {finput.close();loadCompressFile(Fname);return;}
    loadFile(Fname);
    finput.close();
}

Fingerprint::Fingerprint(const string& fgp,const string& name,const bool& numerical) throw(MoleExcept):name(name)
{
    bit_size=0;
    numeric=numerical;
    size_t pos=fgp.find_first_of(":");
    if (pos != string::npos) {loadSVM(fgp);return;}
    pos = fgp.find_first_of("[");
    if (pos != string::npos) {loadCompress(fgp);return;}
    load(fgp,numerical);

}

Fingerprint::~Fingerprint() {}


Fingerprint::Fingerprint(Fingerprint const &FGP)
{
    numeric=FGP.numeric;
    name =FGP.name;
    bit_size=FGP.bit_size;
    bitstring.clear();
    bitstring.reserve(bit_size);
    for (unsigned int i=0; i < bit_size;i++) bitstring.push_back(FGP.at(i));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// MODIFIERS /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
  * \brief Add 1 to the selected bin
  * \fn void Fingerprint::Push_one(const unsigned int _pos)
  * \param _pos : position of the bin
  *
  * If the fingerprint is numeric, it will add one to the bin, otherwise set the bin to 1.
  * \throw MoleExcept 1110101 - _pos is above the size of the fingerprint
  */
void Fingerprint::pushOne(const unsigned int& pos) throw(MoleExcept)
{
    if (pos >= bit_size) throw MoleExcept(1110101,"Fingerprint::pushOne","Given position is above the size of the fingerprint");
    if (numeric) bitstring.at(pos) ++;
    else bitstring.at(pos) =1;

}
/*!
  * \brief Set the given bit on (binary only)
  * \fn void Fingerprint::BitOn   (const unsigned int _pos)
  * \param _pos : position of the bin
  *
  * Set the given bin to 1, only with a binary fingerprint
  * \throw MoleExcept 1110102 - _pos is above the size of the fingerprint
  */
void Fingerprint::bitOn   (const unsigned int& pos)throw(MoleExcept) {
    if (pos >= bit_size) throw MoleExcept(1110201,"Fingerprint::bitOn","Given position is above the size of the fingerprint");
    if(!numeric)bitstring.at(pos)= 1;
}
/*!
  * \brief Set the given bit off (binary only)
  * \fn void Fingerprint::BitOff   (const unsigned int _pos)
  * \param _pos : position of the bin
  *
  * Set the given bin to 0 only with a binary fingerprint
  * \throw MoleExcept 1110103 - _pos is above the size of the fingerprint
  */
void Fingerprint::bitOff  (const unsigned int& pos) throw(MoleExcept){
    if (pos >= bit_size) throw MoleExcept(1110301,"Fingerprint::bitOff","Given position is above the size of the fingerprint");
    if(!numeric)bitstring.at(pos)= 0;
}
/*!
  * \brief Switch the value of the bin (binary only)
  * \fn void Fingerprint::Toggle   (const unsigned int _pos)
  * \param _pos : position of the bin
  *
  * For a binary fingerprint, set the value to 1 if its value is 0 and 0 if the value is 1.
  * \throw MoleExcept 1110104 - _pos is above the size of the fingerprint
  */
void Fingerprint::toggle  (const unsigned int& pos)throw(MoleExcept)
{
    if (pos >= bit_size) throw MoleExcept(1110401,"Fingerprint::toggle","Given position is above the size of the fingerprint");
    if(!numeric) bitstring.at(pos)= (bitstring.at(pos)) ?0:1;
}

/*!
  * \brief Add value to the bin
  * \fn void Fingerprint::Add_data(const unsigned int _pos, const unsigned int _val)
  * \param _pos : position of the bin
  * \param _val : value to add
  *
  *  For a binary fingerprint, set the value to 1 if _val is equal or more than one.<br/>
  * For a numerical fingerprint, add to the bin the given value.
  * \throw MoleExcept 1110105 - _pos is above the size of the fingerprint
  */
void Fingerprint::addData(const unsigned int& pos, const unsigned int& val) throw(MoleExcept){
    if (pos >= bit_size) throw MoleExcept(1110501,"Fingerprint::addData","Given position is above the size of the fingerprint");
    if (numeric) bitstring.at(pos)+=val;
    else if (val >=1)bitstring.at(pos)=1;
}


/*!
  * \brief Toggle all value
  * \fn void Fingerprint::Toggle_All()
  * Same as Toggle but for all bins
  * set the value to 1 if its value is 0 and 0 if the value is above or equal to 1.
  */
void Fingerprint::toggleAll(){
    for (std::vector<unsigned int>::iterator it = bitstring.begin(); it != bitstring.end(); it++ ) (*it) = (*it)? 0:1;
}

/*!
  * \brief Accessing value of a given bin
  * \fn unsigned int& Fingerprint::operator[] (const unsigned int&_pos)
  * \param _pos : position of the bin
  * \warning It can throw an exception if _pos is abouve the length of the fingerprint
  *
    * Return the value of the pos bin
  */
unsigned int& Fingerprint::operator[] (const unsigned int&_pos)throw(MoleExcept) {
    if (_pos >= bit_size) throw MoleExcept(1110601,"Fingerprint::operator[]","Given position is above the size of the fingerprint");
    return bitstring.at(_pos);

}


/*!
  * \brief Convert a numerical fingerprint to binary
  * \fn void Fingerprint::ToBinary()
  *
  * Scan each bin and set them to 1 if there values is equal or more than 1.
  */
void Fingerprint::toBinary()
{
    for (std::vector<unsigned int>::iterator it = bitstring.begin(); it != bitstring.end(); it++ ) (*it) = (*it)? 1:0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// OUTPUT ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
  * \brief Convert the fingerprint to a string
  * \fn std::string Fingerprint::ToString()
  *
  * Return a string of each value
  */
std::string Fingerprint::toString() const
{
    std::ostringstream oss;
    for (vector<unsigned int>::const_iterator it = bitstring.begin(); it != bitstring.end(); it++)
    {
        oss << (*it);
    }
    return oss.str();
}

/*!
  * \brief Convert the fingerprint to a string
  * \fn std::string Fingerprint::ToString()
  *
  * Return a string of each value separated by a space
  */

std::string Fingerprint::toStringbin() const
{
    std::ostringstream oss;
    for (vector<unsigned int>::const_iterator it = bitstring.begin(); it != bitstring.end(); it++)
    {
        oss << (*it) << " " ;
    }
    return oss.str();
}

/*!
  * \brief Convert the fingerprint to a compressed string
  * \fn std::string Fingerprint::CompressString()
  *
  * Return a compressed string. Every value not 0 are outputed. A serie of 0 will be outputed by [ followed by the number of 0.<br/>
  * Example : 0 1 0 0 1 0 1 0 0 0<br/>
  * Return :  [1 1 [2 1 [1 1 [3<br/>
  */
std::string Fingerprint::toCompressString() const
{
    ostringstream oss;
    unsigned int count=0;
    unsigned int inter =0;
    for (vector<unsigned int>::const_iterator it = bitstring.begin(); it != bitstring.end(); it++)
    {
        inter++;
        if ((*it) == 0) count++;
        else {
            if (count !=0){ oss << "[" << count <<"-" << (*it) << "-"; count=0;}
            else oss << (*it)<< "-";
        }

    }
    if (count !=0){ oss << "[" << count <<"-"; count=0;}
    return oss.str();


}


/*!
  * \brief Convert the fingerprint to an SVM string
  * \fn std::string  Fingerprint::SVM_String()
  *
  * Example : 0 1 0 0 1 0 1 0 0 0<br/>
  * Return :  1:0 2:1 3:0 4:0 5:1 6:0 7:1 8:0 9:0 10:0<br/>
  */
string  Fingerprint::toSVMString() const
{
    ostringstream oss;
    unsigned int pos=0;

    for (vector<unsigned int>::const_iterator it = bitstring.begin(); it != bitstring.end(); it++)
    {
        pos++;
        if (pos == bit_size ||(*it) != 0) { oss << pos<<":"<<(*it)<<" ";}
    }
    return oss.str();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// INPUTS ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
  * \brief Transform a compressed string into a fingerprint
  * \fn void Fingerprint::Load_Compress(string str)
  * \param str : compressed string
  */
void Fingerprint::loadCompress(const string & str) throw(MoleExcept)
{
    const unsigned int size = str.length() ;
    if (size == 0) throw MoleExcept(1110701, "Fingerprint::loadCompress","No fingerprint given");
    ostringstream oss;

    char c;
    unsigned int Num_zero=0;
    unsigned int pos=0;
    for (unsigned int Ivar1 = 0; Ivar1 < size; Ivar1++)
    {
        c = str[Ivar1];
        if (c == '[')
        {
            oss.str(""); Ivar1++; c=str[Ivar1]; oss << c;
            while(c != '-') {
                Ivar1++;
                c=str[Ivar1];
                if (c!= '-')oss << c;
            }
            Num_zero = atoi(oss.str().c_str()); pos += Num_zero;
            for (unsigned int List_zero = 0; List_zero< Num_zero; List_zero++) {    bitstring.push_back(0); }
        }
        else
        {
            oss.str("");  oss << c;
            while(c != '-') {
                Ivar1++;
                if (Ivar1 > size) break;
                c=str[Ivar1];
                if (c!= '-')oss << c;
            }
            bitstring.push_back(atoi(oss.str().c_str())); pos++;
        }
    }
    bit_size=pos;
}

/*!
 * \fn void  Fingerprint::load(const string fgp)throw(ICMole::MoleExcept)
  * \brief Load a file to transform the content into a fingerprint
  * \param fgp : fingerprint to analyse
  */
void  Fingerprint::load(const string &fgp, const bool& numerical)throw(ICMole::MoleExcept)
{

    if (!numerical){
       const unsigned int size = fgp.size();
       if (size != bit_size){bit_size = size; bitstring.reserve(size);}
       for (unsigned int Iv = 0; Iv < bit_size; Iv++) bitstring.push_back(0);
       for(std::string::size_type Ivar1 = 0; Ivar1 < fgp.length(); ++Ivar1)
       {
           if (fgp[Ivar1] == '0')bitstring.at(Ivar1) = 0;
           else bitstring.at(Ivar1) = 1;
//           bitstring.at(Ivar1) = stoi(fgp[Ivar1].c_str());
       }
    }else {
        vector<string> values;
        ICMole::tokenStr(fgp, values," ");
        const unsigned int size = values.size();

        if (size ==0) throw MoleExcept(1111201,"Fingerprint::load","Nothing to load");
        if (size != bit_size){bit_size = size; bitstring.reserve(size);}
        for (unsigned int Iv = 0; Iv < bit_size; Iv++) bitstring.push_back(0);

        for (unsigned int Ivar1 =0; Ivar1 < bit_size; Ivar1++)
        {

            bitstring.at(Ivar1) = atoi((values.at(Ivar1)).c_str());

        }
    }
}


void Fingerprint::loadSVM(const string& str) throw(ICMole::MoleExcept)
{
    vector<string> entries;vector<string> posing;
    string bitter;
    ICMole::tokenStr(str,entries," ");

    unsigned int pos;
    bitter = *(entries.begin()+entries.size()-1);
    ICMole::tokenStr(bitter,posing,":");
    bit_size = atoi((*(posing.begin())).c_str());bitstring.reserve(bit_size);

    for (unsigned int Iv = 0; Iv < bit_size; Iv++) bitstring.push_back(0);
    for (vector<string>::iterator its = entries.begin(); its != entries.end(); its++)
    {
        posing.clear();
        bitter = *its;
        ICMole::tokenStr(bitter,posing,":");
        pos = atoi((*(posing.begin())).c_str())-1;
        if (pos >= bit_size)  throw MoleExcept(1110801, "Fingerprint::loadSVM","Given position is above fingerprint length "+bitter);
        bitstring[pos] =atoi((*(posing.begin()+1)).c_str());
    }
}

/*!
 * \fn void  Fingerprint::loadFile(const string _fname)throw(ICMole::MoleExcept)
  * \brief Load a file to transform the content into a fingerprint
  * \param fgp : fingerprint to analyse
  */
void  Fingerprint::loadFile(const string& _fname)throw(ICMole::MoleExcept)
{
    ifstream finput;
    string ligne;
    finput.open(_fname.c_str(),ios::in);
    if (finput.is_open() == false) { throw MoleExcept(1110901,"Fingerprint::LOAD","Cannot open fingerprint file : " + _fname);}

    std::getline(finput,ligne);
    size_t pos=ligne.find_first_of("\t");
    if (pos == string::npos){ throw MoleExcept(1110902,"Fingerprint::LOAD","Cannot read fingerprint file : " + _fname);}
    name = ligne.substr(0,pos);
    ligne = ligne.substr(pos+1);

    finput.close();
    load(ligne);

}

/*!
  * \brief Load a file to transform the content into a fingerprint
  * \param _fname : File to read
  */
void Fingerprint::loadCompressFile(const string& _fname)throw(ICMole::MoleExcept)
{

    ifstream 		finput;
    string ligne;
    finput.open(_fname.c_str(),ios::in);
    if (finput.is_open() == false)  { throw MoleExcept(1111001,"Fingerprint::LOAD_COMPRESS_FILE","Cannot open fingerprint  file : " + _fname);}
    std::getline( finput, ligne );
    size_t pos=ligne.find_first_of("\t");
    if (pos == string::npos){ throw MoleExcept(1111002,"Fingerprint::LOAD","Cannot read fingerprint file : " + _fname);}
    name = ligne.substr(0,pos);
    ligne = ligne.substr(pos+1);
    if (ligne.length() == 0) {throw MoleExcept(1111003,"Fingerprint::LOAD","Nothing to read in " + _fname);}
    finput.close();

    loadCompress(ligne);
}



/*!
   * \brief Load a file to transform the content into a fingerprint
   * \param _fname : File to read
   */
void Fingerprint::loadSVM_File(const string& _fname)throw(ICMole::MoleExcept)
{
    ifstream finput;
    string ligne;

    finput.open(_fname.c_str(),ios::in);
    if (finput.is_open() == false) { throw MoleExcept(1111101,"Fingerprint::loadSVM_File","Cannot open fingerprint svm file : " + _fname);}
    std::getline(finput,ligne);

    size_t posI=ligne.find_first_of("\t");
    if (posI == string::npos){ throw MoleExcept(1111102,"Fingerprint::LOADSVMFile","Cannot read fingerprint file : " + _fname);}
    name = ligne.substr(0,posI);


    ligne = ligne.substr(posI+1);
    if (ligne.length() == 0) {throw MoleExcept(1111103,"Fingerprint::LOAD","Nothing to read in " + _fname);}
    finput.close();

    loadSVM(ligne);
}






