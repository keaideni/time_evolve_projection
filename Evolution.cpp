#include "Evolution.h"
#include <iostream>

using namespace std;


Evolution::Evolution(const MatrixXd& H1, const MatrixXd& H2, const double& t):
_eigenstate(H1),
_eigenstateend(H2)
{
        SelfAdjointEigenSolver<MatrixXd> tepOP(H2);

        //cout<<eigenstate.eigenvalues().size()<<endl;
        //cout<<tepOP.eigenvalues().size()<<endl;

        int n=tepOP.eigenvalues().size();
        VectorXcd timeOP(n);

        for(int i=0; i<n; ++i)timeOP[i]=std::complex<double>(cos(tepOP.eigenvalues()(i)*t),-1*sin(tepOP.eigenvalues()(i)*t));

        MatrixXcd E=timeOP.asDiagonal();
        //cout<<E.adjoint()*E<<endl;
        //cout<<"the real:"<<endl;
        //cout<<tepOP.eigenvalues()<<endl;
        //cout<<"the complex"<<endl;
        //cout<<timeOP<<endl;
        //cout<<"the diagonal:"<<endl;
        //cout<<E<<endl;
        
        //cout<<tepOP.eigenvectors()*tepOP.eigenvectors().inverse()<<endl;;
        _tOP=tepOP.eigenvectors()*E*tepOP.eigenvectors().adjoint();
        //cout<<_tOP<<endl;
        for(int i=0; i<n; ++i)timeOP[i]=std::complex<double>(1,0);
        _t0OP=timeOP.asDiagonal();
        
}