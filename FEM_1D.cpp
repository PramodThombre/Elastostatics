#include <iostream>
#include <eigen/Eigen/Dense>
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


//Linear Shape functions and its derivatives
// double N1(int z) {return (1-z)/2;}
// double N1dash(int z) {return -0.5;}
// double N2(int z) {return (1+z)/2;}
// double N2dash(int z) {return +0.5;}

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
    //if (order==1) 
        //NumNodesPerElem=2;
        float quadPts[2] = {-0.57735,0.57735};   //-1/sqrt(3),1/sqrt(3)  
        int NumQuadPts = sizeof(quadPts)/sizeof(int);
        double WeightquadPts[2] = {1,1};   
    // else if (order==2) 
    //     NumNodesPerElem=3;
    //     float quadPts[3] = {-sqrt(3/5), 0, sqrt(3/5)};   //-1/sqrt(3),1/sqrt(3)  
    //     int NumQuadPts = sizeof(quadPts)/sizeof(int);
    //     double WeightquadPts[3] = {5/9, 8/9, 5/9};   
    //typedef Eigen::Matrix<float, 2, 2> KMat;
    //typedef Eigen::Matrix<float, 2, 1> FMat;
    

    double Le = L/NumOfElems; 
    float zDashx = 2/Le;

    for (int i=0;i<NumOfElems;i++)
    {
        Klocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,NumNodesPerElem);
        Flocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,1);
    }
    
    for (int i=0;i<NumOfElems;i++)
    {
        //std::cout << "NumOfElems=" << NumOfElems << std::endl;;
        for (int j=0;j<NumQuadPts;j++)
        {
            //For Klocal
            //std::cout << "quad=" << quadPts[j] << std::endl;
            Eigen::MatrixXf Klocal(NumNodesPerElem,NumNodesPerElem);
            // for (int m=0;m<NumNodesPerElem;m++)
            // {   
            //     for (int n=0;n<NumNodesPerElem;n++)
            //     {
                    //Klocal(m,n) = N1dash(m,order)*N1dash(n,order);
                    if (order == 1)
                    {
                        Klocal(0,0) = N1dash(1,order)*N1dash(1,order);
                        Klocal(1,0) = N1dash(1,order)* N2dash(1,order);
                        Klocal(0,1) = N2dash(1,order)*N1dash(1,order);
                        Klocal(1,1) = N2dash(1,order)*N2dash(1,order);
                    } 
                    //else{
                    //     Klocal(0,0) = N1dash(1)*N1dash(1);
                    //     Klocal(1,0) = N1dash(1)*N2dash(1);
                    //     Klocal(2,0) = N1dash(1)*N3dash(1);
                    //     Klocal(0,1) = N2dash(1)*N1dash(1);
                    //     Klocal(1,1) = N2dash(1)*N2dash(1);
                    //     Klocal(2,1) = N2dash(1)*N3dash(1);
                    //     Klocal(0,2) = N3dash(1)*N1dash(1);
                    //     Klocal(1,2) = N3dash(1)*N2dash(1);
                    //     Klocal(2,2) = N3dash(1)*N3dash(1);
                    // }
            //     }
            // }
            //Eigen::Matrix<float, 2, 2> b {{N1dash(1)*N1dash(1), N2dash(1)*N1dash(1)},{N1dash(1)* N2dash(1), N2dash(1)*N2dash(1)},};
            //std::cout << "Klocal" << Klocal << std::endl;
            Klocal_arr[i] = Klocal_arr[i] + E*A*Klocal*zDashx*WeightquadPts[j];

            //For Flocal
            Eigen::MatrixXf Flocal(NumNodesPerElem,1);
            Flocal(0,0) = N1(quadPts[j],order);
            Flocal(1,0) = N2(quadPts[j],order);
            //std::cout << "Flocal" << Flocal << std::endl;
            if (bondaryCondition == 4)
            {
                //for f=f_bar*x we wll get x=z*L/6
                //std::cout << "A" << A << " Le/2" << (Le/2) << "WeightquadPts " << WeightquadPts[j] << "quadPts " << quadPts[j] << "NumOfElems " << L/(2*NumOfElems) << std::endl;
                //Area is not required as force is defined per unit length
                float z = quadPts[j];
                float x_mid = (i-1)*Le + (Le/2);
                float x_z = z*(Le/2) + x_mid;
                Flocal_arr[i] = Flocal_arr[i] + Flocal*F*(Le/2)*WeightquadPts[j]*x_z;
            } 
            else
            {
                
                //std::cout << "FlocalTemp " << FlocalTemp << std::endl;
                //std::cout << "Flocal_arr " << Flocal_arr[i] << std::endl;

                Flocal_arr[i] = Flocal_arr[i] + Flocal*F*(Le/2)*WeightquadPts[j];;
            }
                //std::cout << "A" << A << " Le/2" << (Le/2) << "WeightquadPts " << WeightquadPts[j] << "quadPts " << quadPts[j] << "NumOfElems " << L/(2*NumOfElems) << std::endl;
            //std::cout << "Flocal_arr " << Flocal_arr[i] << std::endl;
        }

    }
    for (int j=0;j<NumQuadPts;j++)
    {
        //std::cout << "Klocal_arr " << Klocal_arr[j] << std::endl;
        //std::cout << "Flocal_arr " << Flocal_arr[j] << std::endl;
        std::cout << " " << std::endl;
    }
    // for (int j=0;j<NumNodesTotal;j++)
    // {
    //     std::cout << "Flocal_00 j=" <<j <<" "<< Flocal_arr[j](0,0) << std::endl;
    // }
    //std::cout << Matrix3d(mat.triangularView<Lower>()) << "\n\n";
    int i =0;
    for (int count=0;count<NumOfElems;count++)
    {
        for (int m=0;m<NumNodesPerElem;m++)
        {   
            for (int n=0;n<NumNodesPerElem;n++)
            {
                //std::cout << " m=" << m << " n=" << n << " i=" << i << std::endl;
                
                std::cout << "else" <<  std::endl;
                //std::cout << "Klocal_arr[i](m,n) " << Klocal_arr[i](m,n)  << std::endl;
                KGlobal(i+m,i+n) = KGlobal(i+m,i+n) + Klocal_arr[i](m,n);
                //KGlobal(i+m,i+n) = Klocal_arr[i](m,n);
                //std::cout << "KGlobal m,n = (" << i+m << "," << i+n << ")=" << KGlobal(i+m,i+n) << std::endl;
            }
            FGlobal(i+m,0) = FGlobal(i+m,0) + Flocal_arr[i](m,0);
            //std::cout << "FGlobal m,0 = (" << i+m << ",0)=" << FGlobal(i+m,0) << std::endl;
        }
        i = i + NumNodesPerElem - 1;
    }
    // std::cout << "KGlobal " << std::endl;
    // std::cout << KGlobal << std::endl;
    // std::cout << "The matrix KGlobal is of size "
    //         << KGlobal.rows() << "KGlobal" << KGlobal.cols() << std::endl;

    // std::cout << "FGlobal " << std::endl;
    // std::cout << FGlobal << std::endl;
    // std::cout << "The matrix FGlobal is of size "
    //         << FGlobal.rows() << "FGlobal" << FGlobal.cols() << std::endl;
    return 0;
}

int main()
{
    Eigen::MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
  
    // Eigen::MatrixXf Klocal1(NumNodesPerElem,NumNodesPerElem), Klocal2(NumNodesPerElem,NumNodesPerElem);
    // Eigen::MatrixXf Flocal1(NumNodesPerElem,1), Flocal2(NumNodesPerElem,1);
    // std::cout << "Klocal is of size " << Klocal1.rows() << "x" << Klocal1.cols() << std::endl;
    // std::cout << Klocal1 << std::endl;
    
    int NElems = 3;
    int NumNodesPerElem = 2;//linear
    int order = 1; //Linear=1, Quadratic=2
    int NumNodesTotal = NElems+1;//linear
    if (order==1) 
        int NumNodesTotal = NElems+1;//linear
    else if (order==2) 
        int NumNodesTotal = NElems*2+1;//Quadratic
    int bondaryCondition = 1;
    {
       // BC (i) u(0) = g1, u(L) = g2, f = 0
    std::cout << "Boundary Condition 1" << std::endl;
    Eigen::MatrixXf  KGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf  FGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf  DGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    //For first boundary condition F=0
    float F = 0;    
    GetKandF1DFE(NElems, order, bondaryCondition, F, KGlobal, FGlobal);
    std::cout << "KGlobal" << KGlobal << std::endl;
    std::cout << "FGlobal" << FGlobal << std::endl;
    
    //Apply Dirichlet's BCs
    /* for the following sets of boundary conditions and forcing functions:
    (i) u(0) = g1, u(L) = g2, f = 0
    (ii) u(0) = g1, u(L) = g2, f =  f_tilde
    (iii) u(0) = g1, AEu,x = h at x = L, f =  f_tilde
    (iv) u(0) = g1, u(L) = g2, f =  f_bar*x
    where  f_tilde and f_bar are constant values.*/
    //E = 1011 Pa, A = 10−4 m2,  ̃f = 106 Nm−1,  ̄f = 107 Nm−2, L = 0.1 m, g1 = 0, g2 = 0.001 m, h = 106 N
    
    //BC(i) u(0) = g1 = 0, u(L) = g2 = 0.001 , f = 0
    Eigen::MatrixXf  KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    Eigen::MatrixXf  FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    Eigen::MatrixXf  DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    int numRows = NumNodesTotal-2;
    KGlobal_D = KGlobal.block(NumNodesPerElem-1,NumNodesPerElem-1,numRows,numRows);
    //std::cout << "KGlobal_D" << std::endl;
    //std::cout << KGlobal_D << std::endl;
    FGlobal_D = FGlobal.block(1,0,numRows,1);
    //std::cout << "FGlobal_D" << FGlobal_D << std::endl;
    //std::cout << "KGlobal.block(0,2,1,1)" << KGlobal.block(2,2,1,1) << std::endl;
    FGlobal_D = FGlobal_D - KGlobal.block(1,0,numRows,1)*g1 - KGlobal.block(1,3,numRows,1)*g2;
    //std::cout << "FGlobal_D" << FGlobal_D << std::endl;

    DGlobal_D = KGlobal_D.inverse() * FGlobal_D;
    //std::cout << "DGlobal_D" << DGlobal_D << std::endl;
    DGlobal(0,0) = g1;
    DGlobal.block(1,0,numRows,1) = DGlobal_D;
    DGlobal(NumNodesTotal-1,0) = g2;
    std::cout << "DGlobal" << DGlobal << std::endl;
    }

    {
    std::cout << "Boundary Condition 2" << std::endl;

    //BC(ii) u(0) = g1, u(L) = g2, f =  f_tilde
    Eigen::MatrixXf  KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    Eigen::MatrixXf  FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    Eigen::MatrixXf  DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    
    Eigen::MatrixXf KGlobalBC2 = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf FGlobalBC2 = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf DGlobalBC2 = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    //For first boundary condition F=0
    float F = f_tilde;    
    NElems = 3;
    NumNodesPerElem = 2;//linear
    order = 1; //Linear=1, Quadratic=2
    NumNodesTotal = NElems+1;//linear
    if (order==1) 
        int NumNodesTotal = NElems+1;//linear
    else if (order==2) 
        int NumNodesTotal = NElems*2+1;//Quadratic
    bondaryCondition = 2;
    GetKandF1DFE(NElems, order, bondaryCondition, F, KGlobalBC2, FGlobalBC2);
    std::cout << "KGlobalBC2" << KGlobalBC2 << std::endl;
    std::cout << "FGlobalBC2" << FGlobalBC2 << std::endl;
    KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    int numRows = NumNodesTotal-2;
    KGlobal_D = KGlobalBC2.block(NumNodesPerElem-1,NumNodesPerElem-1,numRows,numRows);
    //std::cout << "KGlobal_D" << std::endl;
    //std::cout << KGlobal_D << std::endl;
    FGlobal_D = FGlobalBC2.block(1,0,numRows,1);
    std::cout << "FGlobal_D" << FGlobal_D << std::endl;
    //std::cout << "KGlobal.block(0,2,1,1)" << KGlobalBC2.block(2,2,1,1) << std::endl;
    FGlobal_D = FGlobal_D - KGlobalBC2.block(1,0,numRows,1)*g1 - KGlobalBC2.block(1,3,numRows,1)*g2;
    std::cout << "FGlobal_D" << FGlobal_D << std::endl;

    DGlobal_D = KGlobal_D.inverse() * FGlobal_D;
    //std::cout << "DGlobal_D" << DGlobal_D << std::endl;
    DGlobalBC2(0,0) = g1;
    DGlobalBC2.block(1,0,numRows,1) = DGlobal_D;
    DGlobalBC2(NumNodesTotal-1,0) = g2;
    std::cout << "DGlobalBC2" << DGlobalBC2 << std::endl;
    }

    {
    std::cout << "Boundary Condition 3" << std::endl;
    //BC(iii) u(0) = g1, AEu,x = h at x = L, f =  f_tilde
    Eigen::MatrixXf  KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    Eigen::MatrixXf  FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    Eigen::MatrixXf  DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    
    Eigen::MatrixXf KGlobalBC3 = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf FGlobalBC3 = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf DGlobalBC3 = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    //For first boundary condition F=0
    float F = f_tilde;    
    NElems = 3;
    NumNodesPerElem = 2;//linear
    order = 1; //Linear=1, Quadratic=2
    NumNodesTotal = NElems+1;//linear
    if (order==1) 
        int NumNodesTotal = NElems+1;//linear
    else if (order==2) 
        int NumNodesTotal = NElems*2+1;//Quadratic

    bondaryCondition = 3;
    GetKandF1DFE(NElems, order, bondaryCondition, F, KGlobalBC3, FGlobalBC3);
    FGlobalBC3(NumNodesTotal-1,0) = FGlobalBC3(NumNodesTotal-1,0) + h;
    std::cout << "KGlobalBC3" << KGlobalBC3 << std::endl;
    std::cout << "FGlobalBC3" << FGlobalBC3 << std::endl;
    KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    int numRows = NumNodesTotal-1;
    KGlobal_D = KGlobalBC3.block(NumNodesPerElem-1,NumNodesPerElem-1,numRows,numRows);
    //std::cout << "KGlobal_D" << std::endl;
    //std::cout << KGlobal_D << std::endl;
    FGlobal_D = FGlobalBC3.block(1,0,numRows,1);
    //std::cout << "FGlobal_D" << FGlobal_D << std::endl;
    //std::cout << "KGlobal.block(0,2,1,1)" << KGlobalBC3.block(2,2,1,1) << std::endl;
    FGlobal_D = FGlobal_D - KGlobalBC3.block(1,0,numRows,1)*g1;
    //std::cout << "FGlobal_D" << FGlobal_D << std::endl;

    DGlobal_D = KGlobal_D.inverse() * FGlobal_D;
    //std::cout << "DGlobal_D" << DGlobal_D << std::endl;
    DGlobalBC3(0,0) = g1;
    DGlobalBC3.block(1,0,numRows,1) = DGlobal_D;
    //DGlobalBC3(NumNodesTotal-1,0) = g2;
    std::cout << "DGlobalBC3" << DGlobalBC3 << std::endl;
    }

    {
    std::cout << "Boundary Condition 4" << std::endl;
    //BC(iv) u(0) = g1, u(L) = g2, f =  f_bar*x

    Eigen::MatrixXf  KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    Eigen::MatrixXf  FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    Eigen::MatrixXf  DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    
    Eigen::MatrixXf KGlobalBC4 = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf FGlobalBC4 = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf DGlobalBC4 = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    //For first boundary condition F=0
    float F = f_tilde;    
    NElems = 3;
    NumNodesPerElem = 2;//linear
    order = 1; //Linear=1, Quadratic=2
    NumNodesTotal = NElems+1;//linear
    if (order==1) 
        int NumNodesTotal = NElems+1;//linear
    else if (order==2) 
        int NumNodesTotal = NElems*2+1;//Quadratic

    bondaryCondition = 4;
    GetKandF1DFE(NElems, order, bondaryCondition, F, KGlobalBC4, FGlobalBC4);
    
    std::cout << "KGlobalBC4" << std::endl;
    std::cout << KGlobalBC4 << std::endl;
    std::cout << "FGlobalBC4" << std::endl;
    std::cout << FGlobalBC4 << std::endl;
    KGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,NumNodesTotal-1);
    FGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    DGlobal_D = Eigen::MatrixXf::Zero(NumNodesTotal-1,1);
    int numRows = NumNodesTotal-2;
    KGlobal_D = KGlobalBC4.block(NumNodesPerElem-1,NumNodesPerElem-1,numRows,numRows);
    //std::cout << "KGlobal_D" << std::endl;
    //std::cout << KGlobal_D << std::endl;
    FGlobal_D = FGlobalBC4.block(1,0,numRows,1);
    //std::cout << "FGlobal_D" << FGlobal_D << std::endl;
    //std::cout << "KGlobal.block(0,2,1,1)" << KGlobalBC4.block(2,2,1,1) << std::endl;
    FGlobal_D = FGlobal_D - KGlobalBC4.block(1,0,numRows,1)*g1 - KGlobalBC4.block(1,3,numRows,1)*g2;
    //std::cout << "FGlobal_D" << FGlobal_D << std::endl;

    DGlobal_D = KGlobal_D.inverse() * FGlobal_D;
    //std::cout << "DGlobal_D" << DGlobal_D << std::endl;
    DGlobalBC4(0,0) = g1;
    DGlobalBC4.block(1,0,numRows,1) = DGlobal_D;
    DGlobalBC4(NumNodesTotal-1,0) = g2;
    std::cout << "DGlobalBC4" << DGlobalBC4 << std::endl;
    }

}
