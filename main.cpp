#include <iostream>
#include "com.h"

int OP::Max;

int main()
{
        ifstream inpara("./data/QNosave.txt");
        if(!inpara)
        {
                cerr<<"the file QNosave.txt doesn't exit!"<<endl;
        }

        JC_Parameter para(inpara);

        inpara.close();

        //cout<<"the Parameter is "<<endl;
        //para.show();


        


        OP::Max=para.ParticleNo();

        com(para);
        
}