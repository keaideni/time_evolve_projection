#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"
#include <fstream>

using namespace std;

void CalParNo(const Hamiltanian& Qubit, const Hamiltanian& Resonator, const Parameter& para, const VectorXcd& wave, std::ofstream& qubitoutfile, std::ofstream& resonatoroutfile);
void CalParNo(const Hamiltanian& Qubit, const Hamiltanian& Resonator, const Parameter& para, const VectorXcd& wave, std::ofstream& qubitoutfile, std::ofstream& resonatoroutfile)
{
        OP temp1, temp2;
        int parti;
        for(int i=0; i<para.LatticeSize(); ++i)
        {

                if(i==0)
                {
                        temp1.time(Qubit.SysCDagR(), Qubit.SysCR());
                        for(int j=1; j<para.LatticeSize(); ++j)
                        {
                                if(j%2==1)temp2.kron(temp1, Resonator.SysEye());
                                else temp2.kron(temp1, Qubit.SysEye());

                                temp1=temp2;
                        }

                        MatrixXd tempMat(temp1.QMat()->at(para.ParticleNo()));

                        
                        //cout<<tempMat.rows()<<"x"<<tempMat.cols()<<endl;

                        qubitoutfile<<(wave.adjoint()*tempMat*wave).real()<<"\t";

                }else
                {
                        temp1=Qubit.SysEye();
                        for(int j=1; j<para.LatticeSize(); ++j)
                        {
                                if(j%2==1)
                                {
                                        if(j==i)
                                        temp2.kron(temp1, Resonator.SysCDagR()*Resonator.SysCR());
                                        else temp2.kron(temp1, Resonator.SysEye());
                                }
                                else 
                                {
                                        if(j==i)temp2.kron(temp1, Qubit.SysCDagR()*Qubit.SysCR());
                                        else temp2.kron(temp1, Qubit.SysEye());
                                }
                                temp1=temp2;
                        }

                        MatrixXd tempMat(temp1.QMat()->at(para.ParticleNo()));

                        
                        if(i%2==0)
                        qubitoutfile<<(wave.adjoint()*tempMat*wave).real()<<"\t";
                        else resonatoroutfile<<(wave.adjoint()*tempMat*wave).real()<<"\t";
                }
        }

        qubitoutfile<<endl;resonatoroutfile<<endl;
}



void calEnt(const Hamiltanian& Qubit, const Hamiltanian& Resonator, const Parameter& para, const VectorXcd& wave, ofstream& entoutfile);
void calEnt(const Hamiltanian& Qubit, const Hamiltanian& Resonator, const Parameter& para, const VectorXcd& wave, ofstream& entoutfile)
{
        OP temp1(Qubit.SysEye()), temp2;
        for(int i=1; i<para.ParticleNo(); ++i)
        {
                if(i%2==1)temp2.kron(temp1, Resonator.SysEye());
                else temp2.kron(temp1, Qubit.SysEye());

                temp1=temp2;
        }
        OP saveOP(temp1);

        if(para.ParticleNo()%2==0)temp1=Qubit.SysEye();
        else temp1=Resonator.SysEye();

        for(int i=para.ParticleNo()+1; i<para.LatticeSize(); ++i)
        {
                if(i%2==1)temp2.kron(temp1, Resonator.SysEye());
                else temp2.kron(temp1, Qubit.SysEye());

                temp1=temp2;
        }

        OP waveOP;
        waveOP.Waveinitial(saveOP, temp1, para.ParticleNo());
        int i(0);
        double ent(0);
        for(auto it=waveOP.QMat()->rbegin(); it!= waveOP.QMat()->rend(); ++it)
        {
                MatrixXcd temp(it->second.rows(),it->second.cols());

                for(int rown=0; rown<it->second.rows(); ++rown)
                {
                        for(int coln=0; coln<it->second.cols(); ++coln)
                        {
                                temp(rown, coln)=wave(i++);
                                if(i>wave.size())
                                {
                                        cerr<<"the wave size is wrong!"<<endl;
                                        exit(1);
                                }
                        }
                }

                JacobiSVD<MatrixXcd> svd(temp, ComputeThinU | ComputeThinV);

                for(int j=0; j<svd.singularValues().size(); ++j)
                {
                        ent-=pow(svd.singularValues()(j),2)*log(pow(svd.singularValues()(j),2));
                }
        }



        entoutfile<<ent<<endl;

}