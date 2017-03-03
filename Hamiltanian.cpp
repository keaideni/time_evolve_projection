#include "Hamiltanian.h"


using namespace std;
Hamiltanian::Hamiltanian(const Hamiltanian&H):
_System(H.System()),
_SysEye(H.SysEye()),
_SysCR(H.SysCR()),
_SysCDagR(H.SysCDagR()),
_SysCL(H.SysCL()),
_SysCDagL(H.SysCDagL())
{}



Hamiltanian::Hamiltanian(const qubittype& inital, const JC_Parameter& para):
_System(SpinZ),
_SysEye(SpinEye),
_SysCR(SpinAnnihilation),
_SysCDagR(SpinCreation),
_SysCL(SpinAnnihilation),
_SysCDagL(SpinCreation)
{
        _System.time(para.Wz()/2);
}


Hamiltanian::Hamiltanian(const resonatortype& initial, const JC_Parameter& para):
_SysEye(para,Eye),
_SysCR(para, Annihilation),
_SysCDagR(para, Creation),
_SysCL(para, Annihilation),
_SysCDagL(para, Creation)
{
        _System.time(_SysCDagR, _SysCR);
        _System.time(para.Wc());
}



void Hamiltanian::show() const
{
        cout<<"the System:"<<endl;
        _System.show();
        cout<<"the SysEye:"<<endl;
        _SysEye.show();
        cout<<"the SysCR:"<<endl;
        _SysCR.show();
        cout<<"the SysCDagR:"<<endl;
        _SysCDagR.show();
        cout<<"the SysCL:"<<endl;
        _SysCL.show();
        cout<<"the SysCDagL:"<<endl;
        _SysCDagL.show();

}


void Hamiltanian::kron(const Hamiltanian& HL, const Hamiltanian& HR, const double& coup)
{
        clear();

        _System.kron(HL._System, HR._SysEye);
        OP systemp;
        systemp.kron(HL.SysEye(), HR.System());
        _System.add(systemp);
        systemp.kron(HL.SysCR(), HR.SysCDagL());
        systemp.time(coup);
        _System.add(systemp);
        systemp.kron(HL.SysCDagR(), HR.SysCL());
        systemp.time(coup);
        _System.add(systemp);




        _SysEye.kron(HL.SysEye(), HR.SysEye());
        _SysCR.kron(HL.SysEye(), HR.SysCR());
        _SysCDagR.kron(HL.SysEye(), HR.SysCDagR());
        _SysCL.kron(HL.SysCL(), HR.SysEye());
        _SysCDagL.kron(HL.SysCDagL(), HR.SysEye());
}



void Hamiltanian::final(const double& coup)
{
        OP systemp;
        systemp.time(_SysCDagL, _SysCR);
        systemp.time(coup);
        _System.add(systemp);

        systemp.time(_SysCL, _SysCDagR);
        systemp.time(coup);
        _System.add(systemp);
}


void Hamiltanian::clear()
{
        _System.clear();
        _SysEye.clear();
        _SysCR.clear();
        _SysCDagR.clear();
        _SysCL.clear();
        _SysCDagL.clear();
}




void Hamiltanian::operator=(const Hamiltanian& H)
{
        _System=H._System;
        _SysCR=H._SysCR;
        _SysEye=H._SysEye;
        _SysCDagR=H._SysCDagR;
        _SysCL=H._SysCL;
        _SysCDagL=H._SysCDagL;
}
