#include <iostream>
#include <eigen/Eigen/Dense>
#include <fstream>
//using Eigen::MatrixXd;

//Inputs initialization
 double E = 1e11 ;//Pa
 double A = 1e-4;// m2
 double f_tilde = 1e6;// Nm−1, 
 double f_bar = 1e7 ;//Nm−2, 
 double L = 0.1;// m,
 double g1 = 0;//, 
 double g2 = 0.001;// m, 
 double h = 1e6 ;//N

//Linear Shape functions and its derivatives
double N1(double z, int order) 
{
    if (order == 1)
        return (1-z)/2;
    else if (order == 2)
        return -(1-z)*z/2;
}
double N1dash(double z, int order) 
{   if (order == 1)
        return -0.5;
    else if (order == 2)
        return -(1-(2*z))/2;
}
double N2(double z, int order) {
    if (order == 1)
        return (1+z)/2;
    else if (order == 2)
        return (1-z)*(1+z);
    }
double N2dash(double z, int order) {
    if (order == 1)
        return +0.5;
    else if (order == 2)
        return -(2*z);
    }
double N3(double z, int order) {
    if (order == 1)
        return (1+z)/2;
    else if (order == 2)
        return (1+z)*z/2;
    }
double N3dash(double z, int order) {
    if (order == 1)
        return +0.5;
    else if (order == 2)
        return (1+(2*z))/2;
    }

int GetKandF1DFE(int NumOfElems, int order, int bondaryCondition, float F, Eigen::MatrixXf  &KGlobal, Eigen::MatrixXf  &FGlobal)
{
    int NumNodesTotal = NumOfElems+1;//linear
    Eigen::MatrixXf Klocal_arr[100] = {}; //initialize to NumOfElems
    Eigen::MatrixXf Flocal_arr[100] = {};
    int NumNodesPerElem = 2;
    if (order==1) 
        NumNodesPerElem=2;  
    else if (order==2) 
         NumNodesPerElem=3;

    float quadPts[2] = {-0.57735,0.57735};   //-1/sqrt(3),1/sqrt(3)  
    int NumQuadPts = sizeof(quadPts)/sizeof(int);
    double WeightquadPts[2] = {1,1};   
    

    double Le = L/NumOfElems; 
    float zDashx = 2/Le;

    for (int i=0;i<NumOfElems;i++)
    {
        Klocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,NumNodesPerElem);
        Flocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,1);
    }
    
    for (int i=0;i<NumOfElems;i++)
    {
        for (int j=0;j<NumQuadPts;j++)
        {
            //For Klocal
            Eigen::MatrixXf Klocal(NumNodesPerElem,NumNodesPerElem);
            if (order == 1)
            {
                Klocal(0,0) = N1dash(1,order)*N1dash(1,order);
                Klocal(1,0) = N1dash(1,order)* N2dash(1,order);
                Klocal(0,1) = N2dash(1,order)*N1dash(1,order);
                Klocal(1,1) = N2dash(1,order)*N2dash(1,order);
            } 
            Klocal_arr[i] = Klocal_arr[i] + E*A*Klocal*zDashx*WeightquadPts[j];

            //For Flocal
            Eigen::MatrixXf Flocal(NumNodesPerElem,1);
            Flocal(0,0) = N1(quadPts[j],order);
            Flocal(1,0) = N2(quadPts[j],order);
            if (bondaryCondition == 4)
            {
                 //Area is not required as force is defined per unit length
                float z = quadPts[j];
                float x_mid = i*Le + (Le/2);
                float x_z = z*(Le/2) + x_mid;
                Flocal_arr[i] = Flocal_arr[i] + Flocal*F*(Le/2)*WeightquadPts[j]*x_z;
            } 
            else
            {
                Flocal_arr[i] = Flocal_arr[i] + Flocal*F*(Le/2)*WeightquadPts[j];
            }
        }
    }
    for (int j=0;j<NumQuadPts;j++)
    {
        std::cout << " " << std::endl;
    }
    
    int i =0;
    for (int count=0;count<NumOfElems;count++)
    {
        for (int m=0;m<NumNodesPerElem;m++)
        {   
            for (int n=0;n<NumNodesPerElem;n++)
            {
                KGlobal(i+m,i+n) = KGlobal(i+m,i+n) + Klocal_arr[i](m,n);
            }
            FGlobal(i+m,0) = FGlobal(i+m,0) + Flocal_arr[i](m,0);
        }
        i = i + NumNodesPerElem - 1;
    }
    return 0;
}

Eigen::MatrixXf GetDGlobal(int NElems, int order, int bondaryCondition, float F)
{
     int NumNodesTotal = NElems+1;//linear
    int NumNodesPerElem = 2;//linear
    if (order==1) 
    {
        NumNodesTotal = NElems+1;//linear
        NumNodesPerElem = 2;
    }
    else if (order==2) 
    {
        NumNodesTotal = NElems*2+1;//Quadratic
        NumNodesPerElem = 3;
    }
    
    Eigen::MatrixXf  KGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf  FGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf  DGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    int numRows = NumNodesTotal-2;//For getting F dirichilets 

    GetKandF1DFE(NElems, order, bondaryCondition, F, KGlobal, FGlobal);
    if (bondaryCondition == 3)
    {
        FGlobal(NumNodesTotal-1,0) = FGlobal(NumNodesTotal-1,0) + h;
        numRows = NumNodesTotal-1;
    } else
    {
        numRows = NumNodesTotal-2;
    }
    std::cout << "KGlobal" << KGlobal << std::endl;
    std::cout << "FGlobal" << FGlobal << std::endl;
    
    Eigen::MatrixXf  KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    Eigen::MatrixXf  FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    Eigen::MatrixXf  DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    
    KGlobal_D = KGlobal.block(NumNodesPerElem-1,NumNodesPerElem-1,numRows,numRows);
    FGlobal_D = FGlobal.block(1,0,numRows,1);
    if (bondaryCondition==3)
    {
        FGlobal_D = FGlobal_D - KGlobal.block(1,0,numRows,1)*g1;
    } else
    {
        FGlobal_D = FGlobal_D - KGlobal.block(1,0,numRows,1)*g1 - KGlobal.block(1,NumNodesTotal-1,numRows,1)*g2;
    }

    DGlobal_D = KGlobal_D.inverse() * FGlobal_D;
    DGlobal(0,0) = g1;
    DGlobal.block(1,0,numRows,1) = DGlobal_D;
    if (bondaryCondition!=3)
        DGlobal(NumNodesTotal-1,0) = g2;
    //std::cout << "DGlobal" << DGlobal << std::endl;
    return DGlobal;
}
int main()
{
     //Apply Dirichlet's BCs
    /* for the following sets of boundary conditions and forcing functions:
    (i) u(0) = g1, u(L) = g2, f = 0
    (ii) u(0) = g1, u(L) = g2, f =  f_tilde
    (iii) u(0) = g1, AEu,x = h at x = L, f =  f_tilde
    (iv) u(0) = g1, u(L) = g2, f =  f_bar*x
    where  f_tilde and f_bar are constant values.*/
    //E = 1011 Pa, A = 10−4 m2,  ̃f = 106 Nm−1,  ̄f = 107 Nm−2, L = 0.1 m, g1 = 0, g2 = 0.001 m, h = 106 N
    
    {
    int NElems = 3;
    int order = 1; //Linear=1, Quadratic=2
    std::cout << "I wrote data to the file\n";
    // ifstream is just like cin:

  
    std::cout << "Boundary Condition 1" << std::endl;
    // BC (i) u(0) = g1, u(L) = g2, f = 0
    //For first boundary condition F=0
    float F = 0;    
    int bondaryCondition = 1;
    Eigen::MatrixXf  DGlobal = GetDGlobal(NElems, order, bondaryCondition, F);
    std::cout << "DGlobal" << DGlobal << std::endl;
    std::ofstream out ("DGlobal_HW2E3BC1.txt");
    out << DGlobal << "\n";
    out.close();

    std::cout << "Boundary Condition 2" << std::endl;
    //BC(ii) u(0) = g1, u(L) = g2, f =  f_tilde
    F = f_tilde;    
    bondaryCondition = 2;
    Eigen::MatrixXf  DGlobalBC2 = GetDGlobal(NElems, order, bondaryCondition, F);

    std::ofstream out1 ("DGlobal_HW2E3BC2.txt");
    out1 << DGlobalBC2 << "\n";
    std::cout << "DGlobalBC2" << DGlobalBC2 << std::endl;
    out1.close();
    //out << "DGlobalBC2" << DGlobalBC2 << "\n";

    std::cout << "Boundary Condition 3" << std::endl;
    //BC(iii) u(0) = g1, AEu,x = h at x = L, f =  f_tilde
    F = f_tilde;    
    bondaryCondition = 3;
    //std::ofstream out ("ElemArrHW2_1Elem3.txt");

    Eigen::MatrixXf  DGlobalBC3 = GetDGlobal(NElems, order, bondaryCondition, F);
    std::cout << "DGlobalBC3" << DGlobalBC3 << std::endl;
    //out << "DGlobalBC3" << DGlobalBC3 << "\n";
    std::ofstream out2 ("DGlobal_HW2E3BC3.txt");
    out2 << DGlobalBC3 << "\n";
    out2.close();

    std::cout << "Boundary Condition 4" << std::endl;
    //BC(iv) u(0) = g1, u(L) = g2, f =  f_bar*x
    F = f_bar;    
    bondaryCondition = 4;
    Eigen::MatrixXf  DGlobalBC4 = GetDGlobal(NElems, order, bondaryCondition, F);
   // std::ofstream out ("ElemArrHW2_1Elem3.txt");
   std::ofstream out3 ("DGlobal_HW2E3BC4.txt");
    out3 << DGlobalBC3 << "\n";
    out3.close();

    std::cout << "DGlobalBC4" << DGlobalBC4 << std::endl;
    //out << "DGlobalBC4" << DGlobalBC4 << "\n";
    }
    int NElems = 100;
    int order = 1; //Linear=1, Quadratic=2
    std::cout << "I wrote data to the file\n";
    // ifstream is just like cin:
  
  {
  
    std::cout << "Boundary Condition 1" << std::endl;
    // BC (i) u(0) = g1, u(L) = g2, f = 0
    //For first boundary condition F=0
    float F = 0;    
    int bondaryCondition = 1;
    Eigen::MatrixXf  DGlobal = GetDGlobal(NElems, order, bondaryCondition, F);
    std::cout << "DGlobal" << DGlobal << std::endl;
//    out << "DGlobal" << DGlobal << "\n";
    std::ofstream out ("DGlobal_HW2E100BC1.txt");
    out << DGlobal << "\n";
    out.close();
    
    std::cout << "Boundary Condition 2" << std::endl;
    //BC(ii) u(0) = g1, u(L) = g2, f =  f_tilde
    F = f_tilde;    
    bondaryCondition = 2;
    Eigen::MatrixXf  DGlobalBC2 = GetDGlobal(NElems, order, bondaryCondition, F);
    std::cout << "DGlobalBC2" << DGlobalBC2 << std::endl;
    //out << "DGlobalBC2" << DGlobalBC2 << "\n";
    std::ofstream out1 ("DGlobal_HW2E100BC2.txt");
    out1 << DGlobalBC2 << "\n";
    out1.close();

    std::cout << "Boundary Condition 3" << std::endl;
    //BC(iii) u(0) = g1, AEu,x = h at x = L, f =  f_tilde
    F = f_tilde;    
    bondaryCondition = 3;
    Eigen::MatrixXf  DGlobalBC3 = GetDGlobal(NElems, order, bondaryCondition, F);
    std::cout << "DGlobalBC3" << DGlobalBC3 << std::endl;
   // out << "DGlobalBC3" << DGlobalBC3 << "\n";
    std::ofstream out2 ("DGlobal_HW2E100BC3.txt");
    out2 << DGlobalBC3 << "\n";
    out2.close();

    std::cout << "Boundary Condition 4" << std::endl;
    //BC(iv) u(0) = g1, u(L) = g2, f =  f_bar*x
    F = f_bar;    
    bondaryCondition = 4;
    Eigen::MatrixXf  DGlobalBC4 = GetDGlobal(NElems, order, bondaryCondition, F);
    std::cout << "DGlobalBC4" << DGlobalBC4 << std::endl;
   // out << "DGlobalBC4" << DGlobalBC4 << "\n";
   std::ofstream out3 ("DGlobal_HW2E100BC4.txt");
    out3 << DGlobalBC4 << "\n";
    out3.close();
  }
    //out.close();
    return 0;
}
