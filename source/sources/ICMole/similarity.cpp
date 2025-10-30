#include <math.h>
#include <iostream>
#include "headers/ICMole/similarity.h"


using namespace ICMole;
using namespace std;

Similarity::Similarity(const Fingerprint& ref, const Fingerprint& comp, const bool& num):numeric(num)
    {
        Ref=&ref;
        Comp=&comp;

    }
Similarity::Similarity(const bool& num):Ref((Fingerprint*)NULL),Comp((Fingerprint*)NULL),numeric(num){}


void Similarity::setRef(const Fingerprint& ref){Ref = &ref;}
void Similarity::setComp(const Fingerprint& comp){ Comp = &comp;}



double Similarity::Hamming() const 
    {
        bool weight=false;
        if (Ref->bit_size != Comp->bit_size) throw MoleExcept(1140101,"Similarity::Hamming","Reference and comparison fingerprint are not the same size");
        if (weights.size() > 0)
         {  if (weights.size() != Ref->bit_size) throw MoleExcept(1140102,"Similarity::Hamming","Missing weights");
            else weight=true;
        }
        const size_t size = Ref->bitstring.size();

        if (numeric)
        {
            int ff1=0,ff2=0;double w=1;
            double vals=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                if (weight) w = weights.at(Ivar1);
                ff1 = Ref->bitstring.at(Ivar1);
                ff2 = Comp->bitstring.at(Ivar1);
                vals+=  w*static_cast<double>(abs(ff1-ff2));
            }
            return vals;
        }
        else
        {

            double sumA=0,sumB=0,sumC=0;

            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {

                const double ff1= static_cast<double>(Ref->bitstring.at(Ivar1));
                const double ff2= static_cast<double>(Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double &w = weights.at(Ivar1);
                    sumA+=w*ff1;sumB+=w*ff2;
                    sumC+=w*ff1*ff2;
                }
                else
                {
                    sumA+=ff1;
                    sumB+=ff2;
                    sumC+=ff1*ff2;
                }

            }
            return (sumA+sumB-2*sumC);
        }
    }

double Similarity::Tanimoto() const
    {
        bool weight=false;
        if (Ref->size() != Comp->size())  throw MoleExcept(1140201,"Similarity::Tanimoto","Reference and comparison fingerprint are not the same size");
        if (weights.size() > 0)
         {  if (weights.size() != Ref->size()) throw MoleExcept(1140202,"Similarity::Tanimoto","Missing weights");
            else weight=true;
        }
        const size_t size = Ref->bitstring.size();

        if (numeric){

            double Asquare=0,Bsquare=0,AB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const double ff1= static_cast<double>(Ref->bitstring.at(Ivar1));
                const double ff2= static_cast<double>(Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    Asquare += w*ff1*ff1;
                    Bsquare += w*ff2*ff2;
                    AB += w*(ff1*ff2);
                }
                else
                {
                    Asquare += ff1*ff1;
                    Bsquare += ff2*ff2;
                    AB += (ff1*ff2);
                }

            }

            if ((Asquare+Bsquare-AB) == 0) return 0;
            return (AB/(Asquare+Bsquare-AB));
        }
        else
        {

            double onlyA=0,onlyB=0, bothAB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const unsigned int& ff1= (Ref->bitstring.at(Ivar1));
                const unsigned int& ff2= (Comp->bitstring.at(Ivar1));


                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    if (ff1 !=0 && ff2 == 0) onlyA+=w;
                    if (ff1 ==0 && ff2 != 0) onlyB+=w;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=w;
                }
                else
                {
                    if (ff1 !=0 && ff2 == 0) onlyA+=1;
                    if (ff1 ==0 && ff2 != 0) onlyB+=1;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=1;
                }

            }

            if (onlyA+onlyB+bothAB==0) return 0;
            else return bothAB/(onlyA+onlyB+bothAB);
        }
    }


double Similarity::RTve() const
    {
        bool weight=false;
        if (Ref->size() != Comp->size())throw MoleExcept(1140301,"Similarity::RTve","Reference and comparison fingerprint are not the same size");
        if (weights.size() > 0)
         {  if (weights.size() != Ref->size()) throw MoleExcept(1140302,"Similarity::RTve","Missing weights");
            else weight=true;
        }
        const size_t size = Ref->bitstring.size();

        if (numeric){

            double Asquare=0,Bsquare=0,AB=0;;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const double ff1= static_cast<double>(Ref->bitstring.at(Ivar1));
                const double ff2= static_cast<double>(Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    Asquare += w*ff1*ff1;
                    Bsquare += w*ff2*ff2;
                    AB += w*(ff1*ff2);
                }
                else
                {
                    Asquare += ff1*ff1;
                    Bsquare += ff2*ff2;
                    AB += (ff1*ff2);
                }

            }

            if ((Asquare+Bsquare-AB) == 0) return 0;
            return (AB/(0.95*Asquare+0.05*Bsquare-AB));
        }
        else
        {
            double onlyA=0,onlyB=0, bothAB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const unsigned int& ff1= (Ref->bitstring.at(Ivar1));
                const unsigned int& ff2= (Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    if (ff1 !=0 && ff2 == 0) onlyA+=w;
                    if (ff1 ==0 && ff2 != 0) onlyB+=w;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=w;
                }
                else
                {
                    if (ff1 !=0 && ff2 == 0) onlyA+=1;
                    if (ff1 ==0 && ff2 != 0) onlyB+=1;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=1;
                }

            }

            if (onlyA+onlyB-bothAB==0) return 0;
            else return bothAB/(0.95*onlyA+0.05*onlyB+bothAB);
        }
    }



double Similarity::FTve() const
    {
        bool weight=false;
        if (Ref->size() != Comp->size()) throw MoleExcept(1140401,"Similarity::FTve","Reference and comparison fingerprint are not the same size");
        if (weights.size() > 0)
         {  if (weights.size() != Ref->size()) throw MoleExcept(1140402,"Similarity::FTve","Missing weights");
            else weight=true;
        }
        const size_t size = Ref->bitstring.size();

        if (numeric){

            double Asquare=0,Bsquare=0,AB=0;;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const double ff1= static_cast<double>(Ref->bitstring.at(Ivar1));
                const double ff2= static_cast<double>(Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    Asquare += w*ff1*ff1;
                    Bsquare += w*ff2*ff2;
                    AB += w*(ff1*ff2);
                }
                else
                {
                    Asquare += ff1*ff1;
                    Bsquare += ff2*ff2;
                    AB += (ff1*ff2);
                }

            }

            if ((Asquare+Bsquare-AB) == 0) return 0;
            return (AB/(0.05*Asquare+0.95*Bsquare-AB));
        }
        else
        {
            double onlyA=0,onlyB=0, bothAB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const unsigned int& ff1= (Ref->bitstring.at(Ivar1));
                const unsigned int& ff2= (Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    if (ff1 !=0 && ff2 == 0) onlyA+=w;
                    if (ff1 ==0 && ff2 != 0) onlyB+=w;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=w;
                }
                else
                {
                    if (ff1 !=0 && ff2 == 0) onlyA+=1;
                    if (ff1 ==0 && ff2 != 0) onlyB+=1;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=1;
                }

            }
            if (onlyA+onlyB-bothAB==0) return 0;
            else return bothAB/(0.05*onlyA+0.95*onlyB+bothAB);
        }
    }

double Similarity::Dice() const
    {
        bool weight=false;
        if (Ref->size() != Comp->size())throw MoleExcept(1140501,"Similarity::Dice","Reference and comparison fingerprint are not the same size");
        if (weights.size() > 0)
         {  if (weights.size() != Ref->size())  throw MoleExcept(1140502,"Similarity::Dice","Missing weights");
            else weight=true;
        }
        const size_t size = Ref->bitstring.size();

        if (numeric)
        {

            double Asquare=0,Bsquare=0,AB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const double ff1= static_cast<double>(Ref->bitstring.at(Ivar1));
                const double ff2= static_cast<double>(Comp->bitstring.at(Ivar1));
                if (weight)
                {
                    const double& w = weights.at(Ivar1);

                AB += w*ff1*ff2;
                Asquare += w*ff1*ff1;
                Bsquare += w*ff2*ff2;
                }
                else
                {
                    AB += ff1*ff2;
                    Asquare += ff1*ff1;
                    Bsquare += ff2*ff2;
                }

            }
            if ((Asquare+Bsquare) == 0 ) return 0;
            return ((2*AB/(Asquare+Bsquare))+1)/2;
        }
        else
        {

            double onlyA=0,onlyB=0, bothAB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const unsigned int& ff1=Ref->bitstring.at(Ivar1);
                const unsigned int& ff2= Comp->bitstring.at(Ivar1);
                if (weight)
                {
                    const double& w = weights.at(Ivar1);
                    if (ff1 !=0 && ff2 == 0) onlyA+=w;
                    if (ff1 ==0 && ff2 != 0) onlyB+=w;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=w;
                }
                else
                {
                    if (ff1 !=0 && ff2 == 0) onlyA+=1;
                    if (ff1 ==0 && ff2 != 0) onlyB+=1;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=1;
                }



            }

            if (onlyA+onlyB==0) return 0;
            else return 2*bothAB/(onlyA+onlyB);
        }
    }


double Similarity::Soergel() const
    {
        bool weight=false;

        if (Ref->size() != Comp->size())throw MoleExcept(1140601,"Similarity::Soergel","Reference and comparison fingerprint are not the same size");
        if (weights.size() > 0)
            {  if (weights.size() != Ref->size()) throw MoleExcept(1140602,"Similarity::Soergel","Missing weights");
               else weight=true;
           }
        const size_t size = Ref->bitstring.size();

        if (numeric)
        {

            double diff=0,maxSum=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const double ff1= static_cast<double>(Ref->bitstring.at(Ivar1));
                const double ff2= static_cast<double>(Comp->bitstring.at(Ivar1));

                if (weight)
                {
                   const double& w = weights.at(Ivar1);
                    diff += w*fabs(ff1-ff2);
                    maxSum += w*((ff1>ff2)?ff1:ff2);
                }
                else
                {
                    diff += fabs(ff1-ff2);
                    maxSum += ((ff1>ff2)?ff1:ff2);
                }


            }
            if (maxSum == 0) return 0;
            return (1-(diff/maxSum));
        }
        else
        {

            double onlyA=0,onlyB=0, bothAB=0;
            for (size_t Ivar1 = 0; Ivar1 < size; Ivar1++)
            {
                const unsigned int& ff1=Ref->bitstring.at(Ivar1);
                const unsigned int& ff2= Comp->bitstring.at(Ivar1);

                if (weight)
                {
                    const double &w = weights.at(Ivar1);
                    if (ff1 !=0 && ff2 == 0) onlyA+=w;
                    if (ff1 ==0 && ff2 != 0) onlyB+=w;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=w;
                }
                else
                {
                    if (ff1 !=0 && ff2 == 0) onlyA+=1;
                    if (ff1 ==0 && ff2 != 0) onlyB+=1;
                    if (ff1 !=0 && ff2 == ff1) bothAB+=1;
                }


            }
            if (onlyA+onlyB-bothAB==0) return 0;
            else return (onlyA+onlyB-2*bothAB)/(onlyA+onlyB-bothAB);
        }
    }

