#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "fingerprint.h"

namespace ICMole
{
/*! \class Similarity
  * \brief Measure Similarity between two Fingerprints
  *
  * It takes two fingerprints, one as reference, and the second as the comparison. Both must have the same size in order to measure a similarity.
  * Moreover, weights can be put on bits.
  */
class Similarity
{
private:
            const Fingerprint*  Ref;       /*!< \brief  Reference fingerprint */
            const Fingerprint*  Comp;      /*!< \brief  Comparison fingerprint */
    std::vector<double>  weights;   /*!< \brief  List of bit weights */
                   bool  numeric;   /*!< \brief  Are the fingerprints numeric */
public:
                         Similarity(const bool& num);  /*!< \brief  Minimal constructor */
                         Similarity(const Fingerprint& ref, const Fingerprint& comp, const bool& numeric=false); /*!< \brief  Standard constructor */
                   void  setRef(const Fingerprint& _f1);  /*!< \brief  Set a new fingerprint as reference */
                   void  setComp(const Fingerprint& _f2); /*!< \brief  Set a new fingerprint as comparison */
    inline         void  setWeight(std::vector<double> t) {weights=t;}  /*!< \brief  Set weights */
                 double  Hamming() const throw(MoleExcept); /*!< \brief  Return the Hamming distance. Goes from infinity (dissimilar) to 0 (similar) for numeric fgp. n (dissimilar) to 0 (similar) for binary fingerprint */
                 double  Tanimoto()const ;/*!< \brief  Return the Tanimoto value. Goes from 0 (dissimilar) to 1 (similar).*/
                 double  RTve() const;    /*!< \brief  Return the Ref Tversky value. Goes from 0 (dissimilar) to 1 (similar). */
                 double  FTve() const;    /*!< \brief  Return the Fit Tversky value. Goes from 0 (dissimilar) to 1 (similar). */
                 double  Dice() const;    /*!< \brief  Return the Dice value. Goes from 0 (dissimilar) to 1 (similar). */
                 double  Soergel() const; /*!< \brief  Return the Soergel value. Goes from 0 (dissimilar) to 1 (similar). */

};
}

#endif // SIMILARITY_H
