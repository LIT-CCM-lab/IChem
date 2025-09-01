#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include "global.h"


namespace ICMole
{
/**
  * \class Fingerprint
  * \brief Fingerprint manager
  * \author DESAPHY Jeremy
  * \version 2.2
  * A fingerprint can be either binary or numerical, as defined in the constructor.<br/>
  * Moreoever it has a fixed length.
  */
class Fingerprint
{
    friend class Similarity;
private:
    std::vector<unsigned int> bitstring;    /*!< array of bins */
            size_t bit_size;     /*!< Length of the fingerprint. Cannot be changed when defined. */
          static unsigned int debug;        /*!< For debug purpose only */
                  std::string name;
                         bool numeric;      /*!< true : numeric fingerprint ; false : binary fingerprint*/

public:
    static inline void runDebug() { Fingerprint::debug=1;}
    ///////////////////
    // CONSTRUCTOR : //
    ///////////////////
    Fingerprint();
    /*!
      * \brief Constructor
      * \fn Fingerprint(const unsigned int size,const bool numeric =false)
      * \param size : Length of the fingerprint
      * \param numeric : Tells whether the fingerprint is binary or numeric -> Default binary
      */
    Fingerprint(const unsigned int& size,const bool& numeric =false);
    Fingerprint(const std::string& name, const bool& numeric =false) throw(MoleExcept);
    Fingerprint(const std::string& fgp,  const std::string& name,const bool& numerical=false) throw(MoleExcept);
    Fingerprint(Fingerprint const &);
    /*!
      * \brief Destructor
      */
    ~Fingerprint();
    ///////////////////
    // MODIFIERS :   //
    ///////////////////
    void  addData (const unsigned int& pos, const unsigned int& val)throw(MoleExcept);
    void  toggle  (const unsigned int& pos)throw(MoleExcept);
    void  bitOn   (const unsigned int& pos)throw(MoleExcept);
    void  bitOff  (const unsigned int& pos)throw(MoleExcept);
    void  pushOne(const unsigned int &pos) throw(MoleExcept);
    unsigned int at(const unsigned int& pos) const {return bitstring.at(pos);}
    void  toggleAll();
    void  toBinary();


    inline  size_t  size() const {return bit_size;}
    unsigned int& operator[] (const unsigned int&b)throw(MoleExcept);
    const std::string& getName() const {return name;}
    void setName(const std::string& N) {name=N;}

    std::string toString() const;
    std::string toStringbin() const;
    std::string toSVMString() const;
    std::string toCompressString() const;

    void load(const std::string &fgp, const bool& numerical=false) throw(MoleExcept);
    void loadCompress    (const std::string& str) throw(MoleExcept);
    void loadSVM(const std::string& str) throw(MoleExcept);
    void loadCompressFile(const std::string& _fname)throw(ICMole::MoleExcept);
    void loadFile        (const std::string& _fname)throw(ICMole::MoleExcept);
    void loadSVM_File    (const std::string& _fname)throw(ICMole::MoleExcept);
    static Fingerprint generateIFP(const std::string& protein_file, const std::string& ligand_file, bool numeric);



    std::vector<unsigned int>::iterator first() {return bitstring.begin();}
    std::vector<unsigned int>::iterator end() {return bitstring.end();}

};
}

#endif
