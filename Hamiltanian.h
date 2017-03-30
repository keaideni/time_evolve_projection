#ifndef HAMILTANIAN_H
#define HAMILTANIAN_H

#include "OP.h"
#include "JC_Parameter.h"


enum qubittype
{
        qubit
};

enum resonatortype
{
        resonator
};


class Hamiltanian
{
private:
        OP _System;
        OP _SysEye;
        OP _SysCR;
        OP _SysCDagR;
        OP _SysCL;
        OP _SysCDagL;
        OP _Edgepro;

public:
        Hamiltanian(){};
        ~Hamiltanian(){};
        Hamiltanian(const Hamiltanian&H);
        Hamiltanian(const qubittype& _qubit, const JC_Parameter& para);//initial for qubit
        Hamiltanian(const resonatortype& _resonator, const JC_Parameter& para);//initial for cavities.

//===============================show====================================
        OP System()const{return _System;};
        OP SysEye()const{return _SysEye;};
        OP SysCR()const{return _SysCR;};
        OP SysCDagR()const{return _SysCDagR;};
        OP SysCL()const{return _SysCL;};
        OP SysCDagL()const{return _SysCDagL;};
        void show()const;


//=============================function==================================
        void kron(const Hamiltanian& HL, const Hamiltanian& HR, const double& coup);
        void final(const double& coup);
        void clear();

        void operator=(const Hamiltanian& H);

};










#endif // HAMILTANIAN_H
