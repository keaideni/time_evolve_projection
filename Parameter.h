#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include <fstream>

#define CLASS_NAME_PRINT(x) #x 

class Parameter
{
private:
        int _LatticeSize;
        int _ParticleNo;
        int _SiteNo;
        int _DeltaQL;
        int _DeltaQR;

        int _D;
        int _SweepNo;
        int _EdgeCondition; //0 for open.

    


public:
        
        double Energy;
        Parameter();
        Parameter(const Parameter&para);
        virtual ~Parameter();
        Parameter(std::ifstream& infile);
        void show() const;
        int LatticeSize() const{return _LatticeSize;};
        int ParticleNo()const{return _ParticleNo;};
        int SiteNo()const {return _SiteNo;};
        int DeltaQL()const {return _DeltaQL;};
        int DeltaQR()const {return _DeltaQR;};
        int D()const {return _D;};
        int SweepNo()const {return _SweepNo;};
        int EdgeCondition()const {return _EdgeCondition;};
       
        
};

#endif