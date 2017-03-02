#ifndef EVOLUTION_H
#define EVOLUTION_H

//#include <Eigen/Dense>
#include "Hamiltanian.h"

//using namespace Eigen;

class Evolution
{
private:
        
public:
        Evolution(){};
        ~Evolution(){};

        SelfAdjointEigenSolver<MatrixXd> _eigenstate;
        SelfAdjointEigenSolver<MatrixXd> _eigenstateend;
        MatrixXcd _tOP;

        MatrixXcd _t0OP;//used for the initial wave, make it complex.




        Evolution(const MatrixXd& H1, const MatrixXd& H2, const double& t);
};











#endif // EVOLUTION_H
