#ifndef JC_PARAMETER_H
#define JC_PARAMETER_H

#include "Parameter.h"

class JC_Parameter:public Parameter
{
private:
        double _Wz;
        double _Wc;
        double _gr;
        double _gl;

public:
        JC_Parameter();
        ~JC_Parameter();

        JC_Parameter(const JC_Parameter&);

        JC_Parameter(std::ifstream& infile);
        void show()const;

        double Wz()const{return _Wz;};
        double Wc()const{return _Wc;};
        double gr()const{return _gr;};
        double gl()const{return _gl;};
        
};






















#endif