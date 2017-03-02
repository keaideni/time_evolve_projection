#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"
#include "Evolution.h"
#include "CalQ.h"

using namespace std;
using namespace Eigen;
void com(const JC_Parameter&para)
{
        
        Hamiltanian Qubit(qubit, para);
        Hamiltanian Resonator(resonator, para);

        Hamiltanian Ham1(Qubit), Ham2;

        

        for(int i=1; i<para.LatticeSize(); ++i)
        {
                if(i%2==0)
                        Ham2.kron(Ham1, Qubit, para.gl());
                else
                        Ham2.kron(Ham1, Resonator, para.gr());
                Ham1=Ham2;
        }
        Ham1.final(para.gl());

        MatrixXd H1=Ham1.System().QMat()->at(para.ParticleNo());

        //cout<<H1.rows()<<"x"<<H1.cols()<<endl;

        Ham1=Qubit;
        for(int i=1; i<para.LatticeSize(); ++i)
        {
                if(i%2==0)
                        Ham2.kron(Ham1, Qubit, para.gr());
                else
                        Ham2.kron(Ham1, Resonator, para.gl());

                Ham1=Ham2;
        }
        Ham1.final(para.gr());

        MatrixXd H2=Ham1.System().QMat()->at(para.ParticleNo());
        //cout<<"first"<<endl;
        //cout<<H2.rows()<<"x"<<H2.cols()<<endl;

        ofstream qubitout("./result/qubit"), resonatorout("./result/resonator");
        ofstream entout("./result/entanglement");
        ofstream difout("./result/difference");
        ofstream energyout("./result/energy");


        Evolution tevo(H1, H2, 1);



        VectorXcd a=tevo._t0OP*tevo._eigenstate.eigenvectors().col(0);

        VectorXcd b;

//=====================================first time======================================================



        //=====================for the difference between waves====================
                
        difout<<abs((a.adjoint()*tevo._eigenstate.eigenvectors().col(0))(0,0))
        <<"\t"<<abs((a.adjoint()*tevo._eigenstateend.eigenvectors().col(0))(0,0))<<endl;
        //=========================================================================

        //===========energy===================
        energyout<<(a.adjoint()*H2*a).real()<<endl;

        //cout<<tevo._tOP.adjoint()*H2*tevo._tOP<<endl<<H2<<endl;


        CalParNo(Qubit, Resonator, para, a, qubitout, resonatorout);

        calEnt(Qubit, Resonator, para, a, entout);
//==================================================================================================

        for(double t=1; t<10; t+=1)
        {
                
                //cout<<tevo._tOP<<endl;

                 b=(tevo._tOP*a);

                 a=b;
//=====================for the difference between waves====================
                
                difout<<abs((a.adjoint()*tevo._eigenstate.eigenvectors().col(0))(0,0))
                <<"\t"<<abs((a.adjoint()*tevo._eigenstateend.eigenvectors().col(0))(0,0))<<endl;
//=========================================================================

                //===========energy===================
                energyout<<abs((a.adjoint()*H2*a)(0,0))<<endl;


                CalParNo(Qubit, Resonator, para, a, qubitout, resonatorout);

                calEnt(Qubit, Resonator, para, a, entout);

        }

        qubitout.close();resonatorout.close();
        entout.close();
        difout.close();
        energyout.close();


        //cout<<a.adjoint()*tevo._eigenstate.eigenvectors().col(1)<<endl;



}