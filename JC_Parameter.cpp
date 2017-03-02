#include "JC_Parameter.h"


JC_Parameter::JC_Parameter(){}

JC_Parameter::~JC_Parameter(){}


JC_Parameter::JC_Parameter(std::ifstream& infile):
Parameter(infile)
{
        std::string str;
        infile >> str >> _Wz
                >> str >> _Wc
                >> str >> _gr
                >> str >> _gl;
}

JC_Parameter::JC_Parameter(const JC_Parameter&para):
Parameter(para),
_Wz(para._Wz),
_Wc(para._Wc),
_gr(para._gr),
_gl(para._gl)
{

}

void JC_Parameter::show()const
{
        Parameter::show();

        std::cout << "Wi= " << _Wz << std::endl
                << "Wc= " << _Wc << std::endl
                << "gr= " << _gr << std::endl
                << "gl= " << _gl << std::endl;
}