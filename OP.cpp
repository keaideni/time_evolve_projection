#include "OP.h"


struct Eigstruct
{
        int q;
        double lamda;
        VectorXd state;
};

std::string itos(int i);
std::string itos(int i)
{
        std::stringstream s;
        s << i;
        return s.str();
}

bool comp(const Eigstruct& a, const Eigstruct& b);
bool comp(const Eigstruct& a, const Eigstruct& b)
{
        return (a.lamda > b.lamda);
}

OP::OP(){}

OP::~OP()
{
        //std::cout<<"the "<<++Class_Num<<"th class destruction!"<<std::endl;
}

OP::OP(const Spin& type)
{



        

        MatrixXd temp(1, 1);
        

        _QDim.insert(std::pair<int, int>(0, 1));
        _QDim.insert(std::pair<int, int>(1, 1));


        switch (type)
        {

                case SpinZ:
                {
                      _RLQ.insert(std::pair<int, int>(0, 0));
                      _RLQ.insert(std::pair<int, int>(1, 1));

                      temp(0, 0) = -1;

                      _QMat.insert(std::pair<int, MatrixXd>(0, temp));


                      temp(0, 0) = 1;
                      _QMat.insert(std::pair<int, MatrixXd>(1, temp));


                      break;
                }
                case SpinCreation:
                {
                      _RLQ.insert(std::pair<int, int>(0, 1));

                      temp(0, 0) = 1;

                      _QMat.insert(std::pair<int, MatrixXd>(0, temp));

                      break;
                }
                case SpinAnnihilation:
                {
                      _RLQ.insert(std::pair<int, int>(1, 0));

                      temp(0, 0) = 1;

                      _QMat.insert(std::pair<int, MatrixXd>(1, temp));


                      break;
                }
                case SpinEye:
                {
                      _RLQ.insert(std::pair<int, int>(0, 0));
                      _RLQ.insert(std::pair<int, int>(1, 1));

                      temp(0, 0) = 1;

                      _QMat.insert(std::pair<int, MatrixXd>(0, temp));
                      _QMat.insert(std::pair<int, MatrixXd>(1, temp));
                      break;
                }
        }

}


OP::OP(const Parameter& para, const OP_modes& type)
{
        int min(para.SiteNo()-para.DeltaQL())
                ,max(para.SiteNo()+para.DeltaQR());

        for(int i=min; i<=max; ++i)_QDim.insert(std::pair<int, int>(i,1));

        switch(type)
        {
                case Eye:
                {
                        for(int i=min; i<=max; ++i)
                        {
                                _RLQ.insert(std::pair<int, int>(i, i));

                                MatrixXd tempm(1,1);
                                tempm << 1;

                                _QMat.insert(std::pair<int, MatrixXd>(i, tempm));
                        }

                        break;
                }
                case Creation:
                {
                        for(int i=min; i<max; ++i)
                        {
                                _RLQ.insert(std::pair<int, int>(i, i+1));
                                MatrixXd tempm(1,1);
                                tempm<<sqrt(i+1);

                                _QMat.insert(std::pair<int, MatrixXd>(i, tempm));

                        }
                        break;
                }
                case Annihilation:
                {
                        for(int i=min+1; i<=max; ++i)
                        {
                                _RLQ.insert(std::pair<int, int>(i, i-1));
                                MatrixXd tempm(1,1);
                                tempm<<sqrt(i);

                                _QMat.insert(std::pair<int, MatrixXd>(i, tempm));
                        }
                        break;
                }
        }
        
}

OP::OP(const OP& a):
_RLQ(a._RLQ),
_QDim(a._QDim),
_QMat(a._QMat)
{

}

OP::OP(const OP& a, const OP& b)
{
    kron(a,b);
    //time(coupling);
}


//=========Kron production====================================
void OP::findDim(const OP& a, const OP& b, std::unordered_map<int, int>& oldDim, std::unordered_map<std::pair<int, int>, int, classcom>& startDim) const
{
        for (int na = 0; na <= OP::Max; ++na)
        {
                auto ita = a._QDim.find(na);
                if(ita == a._QDim.end()) continue;
                for (int nb = 0; nb <= OP::Max; ++nb)
                {
                        auto itb = b._QDim.find(nb);
                        if(itb == b._QDim.end()) continue;
                        int tempQ;
                        tempQ=ita->first + itb->first;

                        bool i = (tempQ <= OP::Max);
                        //bool j = ((a.RLQ.at(ita->first).QN + b.RLQ.at(itb->first).QN) <= Max);

                        if (i)
                        {
                                auto  itc = oldDim.find(tempQ);
                                if (itc == oldDim.end())
                                {
                                        oldDim.insert(std::pair<int, int>(tempQ, ita->second * itb->second));
                                        startDim[std::pair<int, int>(ita->first, itb->first)] = 0;
                                }
                                else
                                {
                                        startDim[std::pair<int, int>(ita->first, itb->first)] = oldDim[itc->first];
                                        oldDim[itc->first] += ita->second * itb->second;
                                }
                        }
                }

        }
}


void OP::kron(const OP& a, const OP&b)
{
        clear();





        //find the dimention of each good quantum number, and label the place to put the kron of two blocks.
        std::unordered_map<std::pair<int, int>, int, classcom> startDim;
        findDim(a, b, _QDim, startDim);



        //first calculate the good quantum number.
        for (auto ita = a._RLQ.begin(); ita != a._RLQ.end(); ita++)
        {
                for (auto itb = b._RLQ.begin(); itb != b._RLQ.end(); itb++)
                {
                        int tempQR;
                        tempQR = ita->first + itb->first;

                        bool i = (tempQR <= Max);
                        bool j = ((a._RLQ.at(ita->first) + b._RLQ.at(itb->first)) <= Max);

                        if (i&&j)
                        {



                                auto itc = _RLQ.find(tempQR);
                                if (itc == _RLQ.end())
                                {
                                        int tempQL;
                                        tempQL = ita->second + itb->second;
                                        _RLQ.insert(std::pair<int, int>(tempQR, tempQL));



                                        int col = _QDim.at(tempQR);
                                        int row = _QDim.at(tempQL);

                                        //creat a zero mat matched the dimention of good quantum number.
                                        MatrixXd tempm(MatrixXd::Zero(row, col));
                                        
                                        _QMat.insert(std::pair<int, MatrixXd>(tempQR, tempm));
                                }

                        }
                }
        }









        //calculate the kron and put it in the place.
        for (auto ita = a._QMat.begin(); ita != a._QMat.end(); ita++)
        {
                for (auto itb = b._QMat.begin(); itb != b._QMat.end(); itb++)
                {




                        int tempQR;
                        tempQR = ita->first + itb->first;

                        bool i = (tempQR <= Max);
                        bool j = ((a._RLQ.at(ita->first) + b._RLQ.at(itb->first) )<= Max);

                        if (i&&j)
                        {


                                

                                


                                int startR = startDim[std::pair<int, int>(ita->first, itb->first)];




                                int startL = startDim[std::pair<int, int>(a._RLQ.at(ita->first), b._RLQ.at(itb->first))];



                                for (int m = 0; m < ita->second.rows(); m++)
                                {
                                        for (int n = 0; n < ita->second.cols(); n++)
                                        {
                                                int a = itb->second.rows();
                                                int b = itb->second.cols();
                                                
                                                //QMat.at(tempQR).block(startL + m*itb->second.rows(), startR + n*itb->second.cols(),a,b)
                                                //      = ita->second(m, n)*itb->second;
                                                for (int i = 0; i < itb->second.rows(); ++i)
                                                {
                                                        for (int j = 0; j < itb->second.cols(); ++j)
                                                        {
                                                                _QMat.at(tempQR)(startL + m*itb->second.rows() + i, startR + n*itb->second.cols() + j)
                                                                                = ita->second(m, n)*itb->second(i,j);
                                                        }
                                                }
                                                
                                        }
                                }



                                
                        }

                }
        }


}


//============================================================



//==========the algorithm of operator=========================
void OP::trans(const OP& a)
{
        _QDim = a._QDim;

        for (auto it = a._RLQ.begin(); it != a._RLQ.end(); it++)
        {
                _RLQ.insert(std::pair<int,int>(it->second, it->first));
                _QMat.insert(std::pair<int, MatrixXd>(it->second, a._QMat.at(it->first).transpose()));
        }
        
}


void OP::add(const OP& a, const OP& b)
{
        if (a._QDim.size() != b._QDim.size())
        {
                std::cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }

        for (auto tempit = a._QDim.begin(); tempit != a._QDim.end(); tempit++)
        {
                if (tempit->second != b._QDim.at(tempit->first))
                {
                        std::cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }
        }

        _QDim = a._QDim;
        _RLQ = a._RLQ;
        _QMat = a._QMat;





        for (auto tempit = b._QMat.begin(); tempit != b._QMat.end(); tempit++)
        {
                auto it = _QMat.find(tempit->first);
                if (it == _QMat.end())
                {
                        _QMat.insert(std::pair<int, MatrixXd>(tempit->first, b._QMat.at(tempit->first)));
                        _RLQ.insert(std::pair<int, int>(tempit->first, b._RLQ.at(tempit->first)));

                }
                else
                {
                        if (a._RLQ.at(tempit->first) != b._RLQ.at(tempit->first))
                        {
                                std::cout << "add0: RQ != LQ, WRONG!" << std::endl;
                                exit(1);
                        }
                        _QMat.at(tempit->first) += b._QMat.at(tempit->first);
                }
        }


}




OP OP::operator+(const OP& a)
{
        OP sum;


        //test the QDim.
        if (_QDim.size() != a._QDim.size())
        {
                std::cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }
        for (auto tempit = _QDim.begin(); tempit != _QDim.end(); tempit++)
        {
                if (tempit->second != a._QDim.at(tempit->first))
                {
                        std::cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }

        }



        //test the RLQ.
        

        sum._QDim = _QDim;
        sum._RLQ = _RLQ;
        sum._QMat = _QMat;




        for (auto tempit = a._QMat.begin(); tempit != a._QMat.end(); tempit++)
        {
                auto it = _QMat.find(tempit->first);
                if (it == _QMat.end())
                {
                        sum._QMat.insert(std::pair<int, MatrixXd>(tempit->first, a._QMat.at(tempit->first)));
                        sum._RLQ.insert(std::pair<int, int>(tempit->first, a._RLQ.at(tempit->first)));
                }

                else
                {
                        if (_RLQ.at(tempit->first) != a._RLQ.at(tempit->first))
                        {
                                std::cerr << "+: RQ != LQ, WRONG!" << std::endl;
                                exit(1);
                        }
                        sum._QMat.at(tempit->first) += a._QMat.at(tempit->first);
                }
        }




        return sum;
}




void OP::add(const OP& a)
{
        if (_QDim.size() != a._QDim.size())
        {
                std::cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }

        for (auto tempit = _QDim.begin(); tempit != _QDim.end(); tempit++)
        {
                if (tempit->second != a._QDim.at(tempit->first))
                {
                        std::cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }
        }
        for (auto tempit = a._QMat.begin(); tempit != a._QMat.end(); tempit++)
        {
                auto it = _QMat.find(tempit->first);
                if (it == _QMat.end())
                {
                        _QMat.insert(std::pair<int, MatrixXd>(tempit->first, a._QMat.at(tempit->first)));
                        _RLQ.insert(std::pair<int, int>(tempit->first, a._RLQ.at(tempit->first)));
                }

                else
                {
                        /*if(RLQ.at(tempit->first).QN != a.RLQ.at(tempit->first).QN)
                        {
                        std::cout<<"add1: RQ != LQ, WRONG!"<<std::endl;
                        exit(1);
                        }*/

                        _QMat.at(tempit->first) += a._QMat.at(tempit->first);

                }
        }

        //return *this;
}




void OP::time(const OP& a, const OP& b)
{
        clear();
        if (a._QDim.size() != b._QDim.size())
        {
                std::cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }
        for (auto tempit = a._QDim.begin(); tempit != a._QDim.end(); tempit++)
        {
                if (tempit->second != b._QDim.at(tempit->first))
                {
                        std::cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }

        }

        _QDim = a._QDim;




        //for the times of two operator, it is a good block systerm. just cross it.
        for (auto ita = b._QMat.begin(); ita != b._QMat.end(); ita++)
        {
                auto it = a._QMat.find(b._RLQ.at(ita->first));
                if (it != a._QMat.end())
                {

                        _RLQ.insert(std::pair<int, int>(ita->first, a._RLQ.at(it->first)));
                        _QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

                }
        }

}



OP OP::operator*(const OP& a)
{
        OP Product;
        if (a._QDim.size() != _QDim.size())
        {
                std::cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }
        for (auto tempit = a._QDim.begin(); tempit != a._QDim.end(); tempit++)
        {
                if (tempit->second != _QDim.at(tempit->first))
                {
                        std::cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }

        }

        Product._QDim = a._QDim;




        //for the times of two operator, it is a good block systerm. just cross it.
        for (auto ita = a._QMat.begin(); ita != a._QMat.end(); ita++)
        {
                auto it = _QMat.find(a._RLQ.at(ita->first));
                if (it != _QMat.end())
                {

                        Product._RLQ.insert(std::pair<int, int>(ita->first, _RLQ.at(it->first)));
                        Product._QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

                }
        }
        return Product;
}




void OP::time(const double& x)
{
        for (auto it = _QMat.begin(); it != _QMat.end(); it++)
        {
                _QMat.at(it->first) = _QMat.at(it->first) * x;
        }

      
}




void OP::time(const OP& a, const double& x)
{
        clear();
        _QDim = a._QDim;
        _RLQ = a._RLQ;
        for (auto it = a._QMat.begin(); it != a._QMat.end(); it++)
        {
            _QMat[it->first] = a._QMat.at(it->first)*x;
        }
}




OP OP::operator*(const double& x)
{
        OP times;


        times._QDim = _QDim;
        times._RLQ = _RLQ;


        //the matrix times the parameter number.
        for (auto tempit = _QMat.begin(); tempit != _QMat.end(); tempit++)
        {
                times._QMat.insert(std::pair<int, MatrixXd>(tempit->first, x * _QMat.at(tempit->first)));
        }

        return times;
}



void OP::time(const double& x, const OP& a)
{
        clear();
        _QDim = a._QDim;
        _RLQ = a._RLQ;
        for (auto it = a._QMat.begin(); it != a._QMat.end(); it++)
        {
                _QMat[it->first] = a._QMat.at(it->first)*x;
        }
}


OP operator*(const double& x, const OP& a)
{
        OP time;

        time._QDim = a._QDim;
        time._RLQ = a._RLQ;
        for (auto it = a._QMat.begin(); it != a._QMat.end(); it++)
        {
                time._QMat[it->first] = a._QMat.at(it->first)*x;
        }

        return time;
}



void OP::operator=(const OP& a)
{
        _RLQ=a._RLQ;
        _QDim=a._QDim;
        _QMat=a._QMat;
}



void OP::save(std::ofstream& outfile)
{
        outfile.precision(30);


        outfile << _QDim.size() << std::endl;
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                outfile << it->first << "       " << it->second << "        ";
        }


        //save the RLQ.
        outfile << _RLQ.size() << std::endl;
        for (auto it = _RLQ.begin(); it != _RLQ.end(); it++)
        {
                outfile << it->first << "          " << it->second << "             ";
        }
        outfile << std::endl;



        //save the QMat.
        for (auto it = _QMat.begin(); it != _QMat.end(); it++)
        {
                outfile << it->first << std::endl;

                outfile<< it->second.rows() << "     " <<it->second.cols()<<std::endl;

                outfile << it->second << std::endl;
                /*for(int i=0; i<it->second.n_rows; i++)
                {
                        for(int j=0; j<it->second.n_cols; j++)
                        {
                                outfile <<std::setprecision(20)<< it->second(i,j)<<std::endl;
                        }
                }*/

        }
}



void OP::read(std::ifstream& infile)
{
        //std::ifstream infile(str);
        /*if(!infile.is_open())
        {
                std::cout<<"the file "<< str << " is not exit " << std::endl;
                exit(1);
        }*/

        //read in the QDim.
        clear();
        int size1;
        infile >> size1;

        int tempQ;
        int tempint;
        for (int i = 0; i < size1; i++)
        {

                infile >> tempQ >> tempint;
                _QDim[tempQ] = tempint;
        }

        //read in the RLQ.
        int size2;
        infile >> size2;
        int tempQ1;
        for (int i = 0; i < size2; i++)
        {
                infile >> tempQ >> tempQ1;
                _RLQ[tempQ] = tempQ1;
        }




        //read in the QMat.
        for (int it = 0; it < _RLQ.size(); ++it)
        {

                int tempQR;

                infile >> tempQR;

                int row, col;

                infile >> row >> col;

                MatrixXd A(row, col);
                for (int i = 0; i < row; i++)
                {
                        for (int j = 0; j < col; j++)
                        {
                                infile >> A(i, j);

                        }
                }

                _QMat[tempQR] = A;
        }
}




//=============================================================

//===================for the OPWave============================
void OP::Waveinitial(const OP& SysEye, const OP& EnvEye, const int& partiNo)
{
        clear();

        for(auto QL=SysEye.QDim()->begin(); QL!=SysEye.QDim()->end(); ++QL)
        {
                for(auto QR=EnvEye.QDim()->begin(); QR!=EnvEye.QDim()->end(); ++QR)
                        {
                            //std::cout<<QR->first<<", "<<QL->first<<std::endl;
                                if(QL->first+QR->first==partiNo)
                                {
                                        _RLQ.insert(std::pair<int, int>(QR->first, QL->first));

                                        MatrixXd temp(MatrixXd::Zero(QL->second, QR->second));

                                        _QMat.insert(std::pair<int, MatrixXd>(QR->first, temp));
                                }
                        }
        }   
}






int OP::ltime(const OP& a, const OP& wave)
{
        clear();
        int flag(0);
        for (auto ita = wave._QMat.begin(); ita != wave._QMat.end(); ita++)
        {
                auto it = a._QMat.find(wave._RLQ.at(ita->first));
                if (it != a._QMat.end())
                {
                        flag++;
                        _RLQ.insert(std::pair<int, int>(ita->first, a._RLQ.at(it->first)));
                        _QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

                }
        }
        return flag;
}



int OP::rtime(const OP& wave, const OP& b)
{
        clear();

        int flag(0);

        for (auto ita = b.QMat()->begin(); ita != b.QMat()->end(); ita++)
        {
                for (auto itQ = wave.RLQ()->begin(); itQ != wave.RLQ()->end(); itQ++)
                {

                        if (ita->first == itQ->first)
                        {
                                flag++;
                                int tempQ(b.RLQ()->at(ita->first));
                                _RLQ[tempQ] = itQ->second;
                                _QMat[tempQ] = (wave._QMat.at(itQ->first)) * ita->second.transpose();

                        }
                }
        }
        return flag;

}

//=============================================================


void OP::clear()
{
        _QDim.clear();
        _QMat.clear();
        _RLQ.clear();
}




void OP::show()const
{
        std::cout << "the QDim: " << std::endl;
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                std::cout << it->first << " => " << it->second << std::endl;
        }
        std::cout << "the RLQ: " << std::endl;
        for (auto it = _RLQ.begin(); it != _RLQ.end(); it++)
        {
                std::cout << it->first << " => " << it->second << std::endl;
        }
        std::cout << "the QMat: " << std::endl;
        for (auto it = _QMat.begin(); it != _QMat.end(); it++)
        {
                std::cout << it->first << " => " << it->second<<std::endl;//.rows()<<"x"<<it->second.cols() << std::endl;
        }
}