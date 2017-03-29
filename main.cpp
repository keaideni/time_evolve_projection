#include <iostream>
#include <fstream>
#include "Mat.h"
#include "mpi.h"

int OP::Max;
using namespace std;

int main(int argc, char* argv[])
{
        ifstream inpara("./data/QNosave.txt");
        if(!inpara)
        {
                cerr<<"the file QNosave.txt doesn't exit!"<<endl;
        }

        JC_Parameter para(inpara);
        OP::Max=para.ParticleNo();

        inpara.close();
        //para.show();

        MPI_Status status;

        int myid, numprocess;

        int groupn(48);
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocess);

        int everygroup(48/numprocess);//cout<<everygroup;
        if(myid==0)
        {
                double gl(0.003), gr(0.027);
                Mat H0(para, gl, gr);
                int Wsize(H0.groundstate().size());

                MPI_Bcast(&Wsize, 1, MPI_INT,
                        0, MPI_COMM_WORLD);
                //cout<<Wsize<<endl;
                
                
                //cout<<"the 0:"<<endl;
                //cout<<H0.groundstate()<<endl;
                MPI_Bcast(&H0.groundstate()(0), Wsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                ofstream overlap("./result/overlap");
                overlap<<gl<<"\t"<<H0.groundstate().adjoint()*H0.groundstate()
                        <<"\t"<<H0.groundstate().adjoint()*H0.MatH()
                                *H0.groundstate()/para.ParticleNo()<<endl;
                for(int i=0; i<everygroup; ++i)
                {
                        gl=0.003+0.0005*(i+1);
                        gr=0.03-gl;
                        Mat H(para, gl, gr);
                        overlap<<gl<<"\t"<<abs(H0.groundstate().adjoint()*H.groundstate())
                                <<"\t"<<(H0.groundstate().adjoint()*H.MatH()
                                *H0.groundstate())/para.ParticleNo()<<endl;

                }
                for(int id=1; id<numprocess; ++id)
                {
                        VectorXd gl1(everygroup);
                        VectorXd overlap1(everygroup);
                        VectorXd energy1(everygroup);
                        
                        MPI_Recv(&gl1(0), everygroup, MPI_DOUBLE, id, id,
                                MPI_COMM_WORLD, &status);
                        
                        //cout<<gl1<<endl;

                        MPI_Recv(&overlap1(0), everygroup, MPI_LONG_DOUBLE, id,
                                id+numprocess, MPI_COMM_WORLD, &status);
                        
                        MPI_Recv(&energy1(0), everygroup, MPI_DOUBLE, id,
                                id+2*numprocess, MPI_COMM_WORLD, &status);
                        //cout<<"from 0:"<<endl;
                        //cout<<overlap1<<endl;

                        for(int i=0; i<everygroup;++i)
                        overlap<<gl1(i)<<"\t"<<overlap1(i)<<"\t"<<energy1(i)<<endl;

                        //delete gl1;
                        //delete energy1;
                        
                }
                overlap.close();

        }
        else
        {
                
                int Wsize;
                MPI_Bcast(&Wsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
                //cout<<Wsize<<endl;
                //double* baseground=new double(Wsize);
                VectorXd baseground(Wsize);
                MPI_Bcast(&baseground(0), Wsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                //cout<<"the 1:"<<endl;
                //cout<<baseground<<endl;
                VectorXd gl(everygroup);
                
                VectorXd overlap(everygroup);
                VectorXd energy(everygroup);


                for(int i=everygroup*myid; i<everygroup*(myid+1);++i)
                {
                        //cout<<i-myid*everygroup<<endl;
                        
                        gl((i-myid*everygroup))=0.003+0.0005*(i+1);
                        double gr=0.03-gl[i-myid*everygroup];
                        Mat H(para, gl[i-myid*everygroup], gr);
                        overlap[i-myid*everygroup]=abs((baseground.adjoint()*H.groundstate()));
                        energy(i-myid*everygroup)=baseground.adjoint()*H.MatH()*baseground;
                        energy(i-myid*everygroup)/=para.ParticleNo();

                        
                        
                }
                //cout<<"from "<<myid<<endl;
                //cout<<overlap<<endl;
                
                

                MPI_Send(&gl(0), everygroup, MPI_DOUBLE, 0, 
                        myid, MPI_COMM_WORLD);
                MPI_Send(&overlap(0), everygroup, MPI_LONG_DOUBLE, 0, 
                        myid+numprocess, MPI_COMM_WORLD);
                MPI_Send(&energy(0), everygroup, MPI_DOUBLE, 0, 
                        myid+2*numprocess, MPI_COMM_WORLD);
                //cout<<myid+numprocess<<endl;
                


                //delete gl;
                //delete energy;
                
                
        }

        


        
        MPI_Finalize();//this is so important. to prove the early exit.

        

        //com(para);
        
}