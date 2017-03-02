#ifndef OP_H
#define OP_H

#include "Parameter.h"
#include <utility>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>



struct classcom
{
        size_t operator()(const std::pair<int, int>& l) const
        {
                return (size_t) ((l.first<<20) ^l.second);
        }
};


enum OP_modes{Eye, Creation, Annihilation};
enum Spin{SpinEye, SpinCreation, SpinAnnihilation, SpinZ};//1/2 spin.

using namespace Eigen;

class OP
{
private:
        std::unordered_map<int, int> _QDim;
        std::unordered_map<int, int> _RLQ;
        std::map<int, MatrixXd> _QMat;

        
public:
        //static int Class_Num;
        static int Max;//used for the maximum quantum number.
        OP();
        ~OP();

        const std::unordered_map<int, int>* QDim()const{return &_QDim;};
        const std::unordered_map<int, int>* RLQ()const{return &_RLQ;};
        const std::map<int, MatrixXd>* QMat()const{return &_QMat;};


//=======initialize the boson operator=========================
        OP(const Spin& type);

        OP(const Parameter& para, const OP_modes& type);

        OP(const OP& a);

        OP(const OP& a, const OP& b);
//=============================================================


//=======Kron production========================================
        void findDim(const OP& a, const OP& b, std::unordered_map<int, int>& oldDim, std::unordered_map<std::pair<int, int>, int, classcom>& startDim) const;
        void kron(const OP& a, const OP&b);
//==============================================================


//=======the algorithm of operator==============================
        void trans(const OP& a);

        void add(const OP& a, const OP& b);

        OP operator+(const OP& a);

        void add(const OP& a);                 //for the operator +=, but don't need to copy anything.

        void time(const OP& a, const OP& b);

        OP operator*(const OP& a);

        void time(const double& x);

        void time(const OP& a, const double& x);

        OP operator*(const double& x);

        void time(const double& x, const OP& a);

        void operator=(const OP& a);

        friend OP operator*(const double& x, const OP& a);
//==============================================================


//==========save and read the operation=========================
        void save(std::ofstream& outfile);
        void read(std::ifstream& infile);
//==============================================================



//==========the OPWave operator=================================
        //OPWave operators have no QDim term and truncate operators do have
        //===========for the OPWave operator==================
        void Waveinitial(const OP& SysEye, const OP& EnvEye, const int& partiNo);
        int ltime(const OP& a, const OP& wave);
        int rtime(const OP& wave, const OP& b);




        void clear();
        void show()const;
};





extern OP operator*(const double& x, const OP& a);


#endif

