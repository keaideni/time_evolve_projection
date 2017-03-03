#ifndef MAT_H
#define MAT_H
#include <Eigen/Core>
#include <SymEigsSolver.h>  // Also includes <MatOp/DenseSymMatProd.h>
#include <iostream>
#include "Hamiltanian.h"

using namespace Eigen;
using namespace Spectra;
class Mat
{
private:
        MatrixXd _MatH;
        VectorXd _groundstate;
public:
        Mat(){};
        ~Mat(){};
        MatrixXd MatH()const{return _MatH;};
        VectorXd groundstate()const{return _groundstate;};



        Mat(const JC_Parameter& para, const double& gl, const double& gr)
        {
                Hamiltanian Qubit(qubit, para);
                Hamiltanian Resonator(resonator, para);

                Hamiltanian Ham1(Qubit), Ham2;
                //Ham1.show();
                

                for(int i=1; i<para.LatticeSize(); ++i)
                {
                        //std::cout<<gl<<std::endl;
                        if(i%2==0)
                                Ham2.kron(Ham1, Qubit, gl);
                        else
                        {
                                        Ham2.kron(Ham1, Resonator, gr);
                                        //Ham2.show();
                        }
                        Ham1=Ham2;

                }
                Ham1.final(gl);


                //int coln(Ham1.System().QMat()->at(para.ParticleNo()).cols());
                //int rown(Ham1.System().QMat()->at(para.ParticleNo()).rows());

                //_MatH.resize(rown, coln);
                //Ham1.System().show();
                _MatH=Ham1.System().QMat()->at(para.ParticleNo());

                DenseSymMatProd<double> opmat(_MatH);

                SymEigsSolver<double, SMALLEST_ALGE, DenseSymMatProd<double>>
                 eigs(&opmat, 1, 4);

                eigs.init();
                eigs.compute();

                _groundstate=eigs.eigenvectors(1);
        };
};








#endif