#include "Parameter.h"

Parameter::Parameter(){}

Parameter::~Parameter(){}

Parameter::Parameter(std::ifstream& infile)
{
        std::string str;
        infile>> str >>_LatticeSize
                >> str >> _ParticleNo
                >> str >> _SiteNo
                >> str >> _DeltaQL
                >> str >> _DeltaQR
                >> str >> _D
                >> str >> _SweepNo
                >> str >> _EdgeCondition
                >> str >> Energy;
}

Parameter::Parameter(const Parameter& para):
_LatticeSize(para._LatticeSize),
_ParticleNo(para._ParticleNo),
_SiteNo(para._SiteNo),
_DeltaQL(para._DeltaQL),
_DeltaQR(para._DeltaQR),

_D(para._D),
_SweepNo(para._SweepNo),
_EdgeCondition(para._EdgeCondition)
{

}


void Parameter::show()const
{
        std::cout << "LatticeSize= " << _LatticeSize << std::endl
                << "PaticleNo= " << _ParticleNo << std::endl
                << "SiteNo= " << _SiteNo << std::endl
                << "DeltaQL= " << _DeltaQL << std::endl
                << "DeltaQR= " << _DeltaQR << std::endl
                << "D= " << _D << std::endl
                << "SweepNo= " << _SweepNo << std::endl
                << "EdgeCondition= " << _EdgeCondition << std::endl;
}

