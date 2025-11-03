#include "headers/ICMole/molecule.h"
using namespace std;
using namespace ICMole;


/*!
 * \brief Chain::Chain
 * \param molecule : parent molecule
 * \param Name : name of the chain
 * \throw 1090101 : No chain name given
 * Standard private constructor that can only be called by a molecule
 */
Chain::Chain(Molecule  &molecule, const std::string& Name) 
    :name(Name),molecule(molecule),inUse(true)
{
    if (Name.length() == 0)
        throw MoleExcept(1090101,
                         "Chain::Chain",
                         "No chain name given");
}






/*!
 * \brief Chain::Chain
 * \param molecule : parent molecule
 * \param Name : name of the chain
 * \throw 1090201 : No parent molecule given
 * \throw 1090202 : No chain name given
 *  Standard private constructor that can only be called by a molecule
 */

Chain::Chain(Molecule  *const molecule,
             const std::string& Name) :
    name(Name),
    molecule(*molecule),inUse(true)
          {
    if (molecule == (Molecule*)NULL)
        throw MoleExcept(1090201,
                         "Chain::Chain",
                         "No parent molecule given");
    if (Name.length() == 0)
        throw MoleExcept(1090202,
                         "Chain::Chain",
                         "No chain name given");

          }






/*!
 * \brief Chain::addResidu
 * \param residu : residu to add
 * \throw 1090301 : No residu given
 * \throw 1090401 : From former residu chain : Residu not given
 * \throw 1090402 : From former residu chain : Residu not found
 * Add the given residu to this chain. This imply to remove this residu
 * from its former chain (if one) and assign to this residu this chain.
 * Finally it add the given residu to the residu list of this chain
 *
 */
void Chain::addResidu(Residu* const residu) 
{
    if (residu== (Residu*)NULL)
        throw MoleExcept(1090301,
                         "Chain::addResidu",
                         "No residu given");

    try
    {
    if (residu->getChain() != (Chain*)NULL && residu->getChain() != this) residu->chain->delResidu(residu);
    }catch (MoleExcept &e)
    {
        e.addTrace("Chain::addResidu ");
        throw;
    }
    residu->chain=this;
    residus.push_back(residu);
}


/*!
 * \brief Chain::delResidu
 * \param residu : residu to remove
 * \throw 1090401 : No residu given
 * \throw 1090402 : Residu not found
 * Remove the given residu from this chain
 * This imply that the given residu is within this chain and not null
 *
 */
void Chain::delResidu(Residu* const residu) 
{
    if (residu == (Residu*)NULL)
        throw MoleExcept(1090401,
                         "Chain::addResidu",
                         "No residu given");

// Searching this residu in the list :
    ItRes it= std::find(residus.begin(),residus.end(),residu);

// Not found : fail
    if (it == residus.end())
        throw MoleExcept(1090402,
                         "Chain::addResidu",
                         "Residu "+residu->getIdentifier()+
                         " not found in Chain "+name+
                         " of molecule "+molecule.getName());
// Otherwise erase :
        residus.erase(it);

}







/*!
 * \brief Chain::toString
 * \return a string describing each residu of this chain
 * getIdentifier is call for each residu of this chain
 */
const string Chain::toString() const
{
    ostringstream oss;
    oss << "###### CHAIN "<< name<<" ##### "<<endl;
    for (ItCRes it = residus.begin(); it != residus.end(); it++)
    {
        oss << (*it)->getIdentifier()<<endl;
    }
    return oss.str();
}







/**
 * @brief Chain::checkUse
 *
 * Check if at least one residu has inUse set to true.
 * When not, the chain inUse value is set to false.
 *
 */
void Chain::checkUse()
{

    for (ItCRes itR= residus.begin(); itR != residus.end(); itR++)
    {
        if ((*itR)->isUsed()) {
            inUse=true;

        }
    }
    inUse=false;

}



/**
 * @brief Chain::setUse
 * @param use: new inUse value
 * @param applyto : when true, apply use value to residu atom and bonds within it
 *
 * Set the new use to the chain and eventually to residu, atom and bonds.
 * When use is set to false, grids, interaction and other tools will ignore it.
 *
 */
void Chain::setUse(const bool& use,const bool& applyto)
{
    inUse=use;
    if (!applyto) return;
    for (ItRes itR = residus.begin();
               itR!= residus.end();
               itR++)
    {
        Residu &res = **itR;

        res.setUse(use,true,true);
    }
}



Residu& Chain::getResidu(const size_t& pos) const 
{
    if (pos >= residus.size()) throw MoleExcept(1071401,"Chain::getResidu","Position is above the number of residus within the chain");
    else return *residus.at(pos);

}

















