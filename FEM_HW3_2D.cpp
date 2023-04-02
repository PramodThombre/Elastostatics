#include <iostream>
#include <eigen/Eigen/Dense>
#include <fstream>
//using Eigen::MatrixXd;

using std::ostringstream; // <-- Added
using Eigen::EigenBase;   // <-- Added

template <typename Derived>
std::string get_shape(const EigenBase<Derived>& x)
{
    std::ostringstream oss;
    oss  << "(" << x.rows() << ", " << x.cols() << ")";
    return oss.str();
}

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

Eigen::MatrixXf GetShapeFn(double z, double eta, int order)
{
    Eigen::MatrixXf N_ShapeFn = Eigen::MatrixXf::Zero(1,4);
    N_ShapeFn(0,0) = (1+z)*(1+eta)/2;
    N_ShapeFn(0,1) = (1-z)*(1+eta)/2;
    N_ShapeFn(0,2) = (1-z)*(1-eta)/2;
    N_ShapeFn(0,3) = (1+z)*(1-eta)/2;
    return N_ShapeFn;
}
Eigen::MatrixXf GetShapeFnGradientwZ(double z, double eta, int order)
{
    Eigen::MatrixXf N_ShapeFnGradient = Eigen::MatrixXf::Zero(1,4);
    N_ShapeFnGradient(0,0) = (1+eta)/2;
    N_ShapeFnGradient(0,1) = -(1+eta)/2;
    N_ShapeFnGradient(0,2) = -(1-eta)/2;
    N_ShapeFnGradient(0,3) = (1-eta)/2;
    return N_ShapeFnGradient;
}
Eigen::MatrixXf GetShapeFnGradientwEta(double z, double eta, int order)
{
    Eigen::MatrixXf N_ShapeFnGradient = Eigen::MatrixXf::Zero(1,4);
    N_ShapeFnGradient(0,0) = (1+z)/2;
    N_ShapeFnGradient(0,1) = (1-z)/2;
    N_ShapeFnGradient(0,2) = -(1-z)/2;
    N_ShapeFnGradient(0,3) = -(1+z)/2;
    return N_ShapeFnGradient;
}


Eigen::MatrixXd GetElemArr(int NHorizontalElems, int NVerticalElems, int NNodesPerElem)
{
    //Elem_LocaltoGlobalNodeArr(ElemId,LocalNodeId) = GlobalNodeId;
    Eigen::MatrixXd Elem_LocaltoGlobalNodeArr = Eigen::MatrixXd::Zero(NHorizontalElems*NVerticalElems,NNodesPerElem);
    int count = 0;
    for(int i=0;i<NVerticalElems;i++)
    {
        //std::cout << "i," << i << std::endl;
        for(int j=0;j<NHorizontalElems;j++)
        {
//            std::cout << "j," << j << std::endl;
                int nElemId = i*NHorizontalElems+j;
                int nGNodeId1 = nElemId+i;
                int nGNodeId2 = nElemId+1+i;
                int nGNodeId3 = nElemId+NHorizontalElems+i+2;
                int nGNodeId4 = nElemId+NHorizontalElems+i+1;
                int k=0;
                //std::cout << "nElemId," << nElemId << std::endl;
                Elem_LocaltoGlobalNodeArr(nElemId,k) = nGNodeId1;
                Elem_LocaltoGlobalNodeArr(nElemId,k+1) = nGNodeId2;
                Elem_LocaltoGlobalNodeArr(nElemId,k+2) = nGNodeId3;
                Elem_LocaltoGlobalNodeArr(nElemId,k+3) = nGNodeId4;
                // std::cout << "nElemId," << nElemId << " nGNodeId1," << nGNodeId1 << std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId2," << nGNodeId2 << std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId3," << nGNodeId3<< std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId4," << nGNodeId4 << std::endl;
                count++;
                //std::cout << "count," << count << std::endl;
        }
    }
    return Elem_LocaltoGlobalNodeArr;
}


Eigen::MatrixXd GetTriaElemArr(int NHorizontalElems, int NVerticalElems, int NNodesPerElem)
{
    //Elem_LocaltoGlobalNodeArr(ElemId,LocalNodeId) = GlobalNodeId;
    Eigen::MatrixXd Elem_LocaltoGlobalNodeArr = Eigen::MatrixXd::Zero(NHorizontalElems*NVerticalElems,NNodesPerElem);
    int count = 0;
    int nElemId = 0;

    for(int i=0;i<NVerticalElems/2;i++)
    {
        //std::cout << "i," << i << std::endl;
        for(int j=0;j<NHorizontalElems/2;j++)
        {
            //std::cout << "j," << j << std::endl;
            
                int nElemIdMultiplier = i*(NHorizontalElems/2)+j;
                int nGNodeId1 = nElemIdMultiplier+i;
                int nGNodeId2 = nElemIdMultiplier+1+i;
                int nGNodeId3 = nElemIdMultiplier+(NHorizontalElems/2)+i+1;
                int k=0;
                //std::cout << "nElemId," << nElemId << std::endl;
                Elem_LocaltoGlobalNodeArr(nElemId,k) = nGNodeId1;
                Elem_LocaltoGlobalNodeArr(nElemId,k+1) = nGNodeId2;
                Elem_LocaltoGlobalNodeArr(nElemId,k+2) = nGNodeId3;
                // std::cout << "nElemId," << nElemId << " nGNodeId1," << nGNodeId1 << std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId2," << nGNodeId2 << std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId3," << nGNodeId3<< std::endl;
            
                nElemId = nElemId+1;
                nElemIdMultiplier = i*(NHorizontalElems/2)+j;
                int nGNodeId11 = nGNodeId2;
                int nGNodeId22= nGNodeId3+1;
                int nGNodeId33 = nGNodeId3;
                k=0;
                //std::cout << "nElemId," << nElemId << std::endl;
                Elem_LocaltoGlobalNodeArr(nElemId,k) = nGNodeId11;
                Elem_LocaltoGlobalNodeArr(nElemId,k+1) = nGNodeId22;
                Elem_LocaltoGlobalNodeArr(nElemId,k+2) = nGNodeId33;
                // std::cout << "nElemId," << nElemId << " nGNodeId11," << nGNodeId11 << std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId22," << nGNodeId22 << std::endl;
                // std::cout << "nElemId," << nElemId << " nGNodeId33," << nGNodeId33<< std::endl;
                nElemId = nElemId+1;
                count++;
                //std::cout << "count," << count << std::endl;
        }
    }
    return Elem_LocaltoGlobalNodeArr;
}

Eigen::MatrixXf GetGNodeCordsArr(int NHorizontalElems, int NVerticalElems, int NNodesPerElem, int dim, float xMin, float yMin,float xMax, float yMax)
{
    float xLength = xMax - xMin;
    float yLength = yMax - yMin;
    float xLe = xLength/NHorizontalElems; //length of each elem in horizontal direction
    float yLe = yLength/NVerticalElems; //length of each elem in vertical direction
    int nNodesHorizontal = NHorizontalElems+1;
    int nNodesVertical = NVerticalElems+1;
    //Elem_LocaltoGlobalNodeArr(ElemId,LocalNodeId) = GlobalNodeId;
    Eigen::MatrixXf Elem_GlobalNodeCordsArr = Eigen::MatrixXf::Zero(nNodesHorizontal*nNodesVertical,dim);
    int count = 0;
    for(int i=0;i<NVerticalElems;i++)
    {
        //std::cout << "i," << i << std::endl;
        for(int j=0;j<NHorizontalElems;j++)
        {
            //std::cout << "j," << j << std::endl;
                int nElemId = i*NHorizontalElems+j;
                int nGNodeId1 = nElemId+i;
                int nGNodeId2 = nElemId+1+i;
                int nGNodeId3 = nElemId+NHorizontalElems+i+2;
                int nGNodeId4 = nElemId+NHorizontalElems+i+1;
                int k=0;
                if(nElemId==399)
                {
                    int debug=1;
                }
                //std::cout << "nElemId," << nElemId << std::endl;
                Elem_GlobalNodeCordsArr(nGNodeId1,k) = xLe*j;//x coord
                Elem_GlobalNodeCordsArr(nGNodeId1,k+1) = yLe*i;//y coord

                Elem_GlobalNodeCordsArr(nGNodeId2,k) = xLe*(j+1);//x coord
                Elem_GlobalNodeCordsArr(nGNodeId2,k+1) = yLe*i;//y coord

                Elem_GlobalNodeCordsArr(nGNodeId3,k) = xLe*(j+1);//x coord
                Elem_GlobalNodeCordsArr(nGNodeId3,k+1) = yLe*(i+1);//y coord

                Elem_GlobalNodeCordsArr(nGNodeId4,k) = xLe*j;//x coord
                Elem_GlobalNodeCordsArr(nGNodeId4,k+1) = yLe*(i+1);//y coord

                // std::cout << "nGNodeId1," << nGNodeId1 << " xCoordGNodeId1," << xLe*j << std::endl;
                // std::cout << "nGNodeId1," << nGNodeId1 << " yCoordGNodeId1," << yLe*i << std::endl;
                // std::cout << "nGNodeId2," << nGNodeId2 << " xCoordGNodeId2," << xLe*(j+1) << std::endl;
                // std::cout << "nGNodeId2," << nGNodeId2 << " yCoordGNodeId2," << yLe*i << std::endl;
                // std::cout << "nGNodeId3," << nGNodeId3 << " xCoordGNodeId3," << xLe*(j+1) << std::endl;
                // std::cout << "nGNodeId3," << nGNodeId3 << " yCoordGNodeId3," << yLe*(i+1) << std::endl;
                // std::cout << "nGNodeId4," << nGNodeId4 << " xCoordGNodeId4," << xLe*j << std::endl;
                // std::cout << "nGNodeId4," << nGNodeId4 << " yCoordGNodeId4," << yLe*(i+1) << std::endl;
                count++;
                //std::cout << "count," << count << std::endl;
        }
    }
    return Elem_GlobalNodeCordsArr;
}

Eigen::MatrixXf GetTriaGNodeCordsArr(int NHorizontalElems, int NVerticalElems, int NNodesPerElem, int dim, float xMin, float yMin,float xMax, float yMax)
{
    float xLength = xMax - xMin;
    float yLength = yMax - yMin;
    float xLe = xLength/(NHorizontalElems/2); //length of each elem in horizontal direction
    float yLe = yLength/(NVerticalElems/2); //length of each elem in vertical direction
    //float sLe = sqrt(xLe^2+yLe^2); // length of each elem in slant direction
    int nNodesHorizontal = (NHorizontalElems/2)+1;
    int nNodesVertical = (NVerticalElems/2)+1;
    //Elem_LocaltoGlobalNodeArr(ElemId,LocalNodeId) = GlobalNodeId;
    Eigen::MatrixXf Elem_GlobalNodeCordsArr = Eigen::MatrixXf::Zero(nNodesHorizontal*nNodesVertical,dim);
    int count = 0;
    int nElemId = 0;
    int divider = (NHorizontalElems/2)+1;
    for(int i=0;i<NVerticalElems/2;i++)
    {
        //std::cout << "i," << i << std::endl;
        for(int j=0;j<NHorizontalElems/2;j++)
        {
            //std::cout << "j," << j << std::endl;
                
                int nElemIdMultiplier = i*(NHorizontalElems/2)+j;
                int nGNodeId1 = nElemIdMultiplier+i;
                int nGNodeId2 = nElemIdMultiplier+1+i;
                int nGNodeId3 = nElemIdMultiplier+(NHorizontalElems/2)+i+1;
                int k=0;
                if(nElemId==399)
                {
                    int debug=1;
                }
                //std::cout << "nElemId," << nElemId << std::endl;
                if(nGNodeId1%divider==0)
                { Elem_GlobalNodeCordsArr(nGNodeId1,k) = xLe*0;}
                else {Elem_GlobalNodeCordsArr(nGNodeId1,k) = xLe*j;}
                Elem_GlobalNodeCordsArr(nGNodeId1,k+1) = yLe*i;//y coord

                Elem_GlobalNodeCordsArr(nGNodeId2,k) = xLe*(j+1);//x coord
                Elem_GlobalNodeCordsArr(nGNodeId2,k+1) = yLe*i;//y coord

                if(nGNodeId3%divider==0)
                { Elem_GlobalNodeCordsArr(nGNodeId3,k) = xLe*0;}//x coord}
                else {Elem_GlobalNodeCordsArr(nGNodeId3,k) = xLe*j;}//x coord}
                Elem_GlobalNodeCordsArr(nGNodeId3,k+1) = yLe*(i+1);//y coord

                std::cout << "nGNodeId1," << nGNodeId1 << " xCoordGNodeId1," << xLe*j << std::endl;
                std::cout << "nGNodeId1," << nGNodeId1 << " yCoordGNodeId1," << yLe*i << std::endl;
                std::cout << "nGNodeId2," << nGNodeId2 << " xCoordGNodeId2," << xLe*(j+1) << std::endl;
                std::cout << "nGNodeId2," << nGNodeId2 << " yCoordGNodeId2," << yLe*i << std::endl;
                std::cout << "nGNodeId3," << nGNodeId3 << " xCoordGNodeId3," << xLe*(j+1) << std::endl;
                std::cout << "nGNodeId3," << nGNodeId3 << " yCoordGNodeId3," << yLe*(i+1) << std::endl;
               
                nElemId = nElemId+1;
                int nGNodeId11 = nGNodeId2;
                int nGNodeId22= nGNodeId3+1;
                int nGNodeId33 = nGNodeId3;
                k=0;

                if(nGNodeId11%divider==0)
                { Elem_GlobalNodeCordsArr(nGNodeId11,k) = xLe*(j+1);}//x coord}
                else {Elem_GlobalNodeCordsArr(nGNodeId11,k) = xLe*j;}//x coord}
                Elem_GlobalNodeCordsArr(nGNodeId11,k+1) = yLe*i;//y coord

                Elem_GlobalNodeCordsArr(nGNodeId22,k) = xLe*(j+1);//x coord
                Elem_GlobalNodeCordsArr(nGNodeId22,k+1) = yLe*(i+1);//y coord

                if(nGNodeId33%divider==0)
                { Elem_GlobalNodeCordsArr(nGNodeId33,k) = xLe*0;}//x coord}
                else {Elem_GlobalNodeCordsArr(nGNodeId33,k) = xLe*j;}//x coord}
                Elem_GlobalNodeCordsArr(nGNodeId33,k+1) = yLe*(i+1);//y coord

                std::cout << "nGNodeId11," << nGNodeId1 << " xCoordGNodeId1," << xLe*j << std::endl;
                std::cout << "nGNodeId11," << nGNodeId1 << " yCoordGNodeId1," << yLe*i << std::endl;
                std::cout << "nGNodeId22," << nGNodeId2 << " xCoordGNodeId2," << xLe*(j+1) << std::endl;
                std::cout << "nGNodeId22," << nGNodeId2 << " yCoordGNodeId2," << yLe*i << std::endl;
                std::cout << "nGNodeId33," << nGNodeId3 << " xCoordGNodeId3," << xLe*(j+1) << std::endl;
                std::cout << "nGNodeId33," << nGNodeId3 << " yCoordGNodeId3," << yLe*(i+1) << std::endl;
               
                nElemId = nElemId+1;
            
                count++;
                //std::cout << "count," << count << std::endl;
        }
    }
    return Elem_GlobalNodeCordsArr;
}

int GetQuadKMandF2DFE(int NumOfElems, Eigen::MatrixXf GNodeCordsArr, Eigen::MatrixXd  ElemArrQuad, int order, int bondaryCondition, float F, Eigen::MatrixXf  &KGlobal, Eigen::MatrixXf  &MGlobal, Eigen::MatrixXf  &FGlobal)
{
    
    int NElemsH =20;
    int NElemsV =20;
    int NPerElems =4;
    int dim = 2;
    float xMin = 0.0;
    float yMin = 0.0;
    float xMax = 1.0;
    float yMax = 1.0;
    double kbar = 385;
    double rho = 6000;
    float rhoC = 3.8151e6; 
    Eigen::MatrixXf Kbar_Mat = Eigen::MatrixXf::Identity(dim,dim);
    Kbar_Mat = kbar*Kbar_Mat;

   
    int NumNodesTotal=(NElemsH+1)*(NElemsV+1);
    Eigen::MatrixXf Klocal_arr[400] = {}; //initialize to NumOfElems
    Eigen::MatrixXf Mlocal_arr[400] = {}; //initialize to NumOfElems
    Eigen::MatrixXf Flocal_arr[400] = {};
    int NumNodesPerElem = 4;
    if (order==1) 
        NumNodesPerElem=4;  
    else if (order==2) 
         NumNodesPerElem=8;

    float quadPts[2] = {-0.57735,0.57735};   //-1/sqrt(3),1/sqrt(3)  
    int NumQuadPts = sizeof(quadPts)/sizeof(int);
    double WeightquadPts[2] = {1,1};   
    
    double L = 1.0;
    double Le = L/NumOfElems; 
    float zDashx = 2/Le;

    for (int i=0;i<NumOfElems;i++)
    {
        Klocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,NumNodesPerElem);
        Mlocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,NumNodesPerElem);
        Flocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,1);
    }
    
   std::cout << "NumOfElems " << NumOfElems << std::endl;    

    for (int i=0;i<NumOfElems;i++)
    {
        for (int j=0;j<NumQuadPts;j++)
        {
            for (int k=0;k<NumQuadPts;k++)
            {
                //For Klocal
                Eigen::MatrixXf Klocal(NumNodesPerElem,NumNodesPerElem);
                if (order == 1)
                {
                    Eigen::MatrixXf N = GetShapeFn(quadPts[j],quadPts[k],1);
                    Eigen::MatrixXf NwZ = GetShapeFnGradientwZ(quadPts[j],quadPts[k],1);
                    Eigen::MatrixXf NwEta = GetShapeFnGradientwEta(quadPts[j],quadPts[k],1);

                    int nodeid1 = ElemArrQuad(i,0);
                    int nodeid2 = ElemArrQuad(i,1);
                    int nodeid3 = ElemArrQuad(i,2);
                    int nodeid4 = ElemArrQuad(i,3);

                    Eigen::MatrixXf nodeIds(1, 4);
                    nodeIds << nodeid1, nodeid2, nodeid3, nodeid4; //comma initialization
                    //std::cout << "nodeIds" << nodeIds << std::endl;

                    Eigen::MatrixXf xCoords(4, 1);
                    xCoords << GNodeCordsArr(nodeid1,0), GNodeCordsArr(nodeid2,0), GNodeCordsArr(nodeid3,0), GNodeCordsArr(nodeid4,0); //comma initialization
                    //std::cout << "xCoords" << xCoords << std::endl;
            
                     Eigen::MatrixXf yCoords(4, 1);
                    yCoords << GNodeCordsArr(nodeid1,1), GNodeCordsArr(nodeid2,1), GNodeCordsArr(nodeid3,1), GNodeCordsArr(nodeid4,1); //comma initialization
                    //std::cout << "yCoords" << yCoords << std::endl;

                    Eigen::MatrixXf Jacobian = Eigen::MatrixXf::Zero(dim,dim);
                    Jacobian << NwZ*xCoords, NwEta*xCoords, NwZ*yCoords, NwEta*yCoords;
                    //std::cout << "Jacobian" << Jacobian << std::endl;
                    
                    Eigen::MatrixXf invJacobian = Eigen::MatrixXf::Zero(dim,dim);
                    invJacobian = Jacobian.inverse();
                    //std::cout << "invJacobian" << invJacobian << std::endl;

                    Eigen::MatrixXf NdashZEta = Eigen::MatrixXf::Zero(NumNodesPerElem,dim);
                    NdashZEta.col(0) << NwZ(0,0),NwZ(0,1),NwZ(0,2),NwZ(0,3);
                    NdashZEta.col(1) << NwEta(0,0),NwEta(0,1),NwEta(0,2),NwEta(0,3);
                    //std::cout << "NdashZEta" << NdashZEta << std::endl;

                    Eigen::MatrixXf NdashX = Eigen::MatrixXf::Zero(NumNodesPerElem,dim);
                    NdashX = NdashZEta*invJacobian;
                    // std::cout << "NdashX" << NdashX << std::endl;
                    // std::cout << "Kbar_Mat" << Kbar_Mat << std::endl;

                    Klocal_arr[i] = Klocal_arr[i] + NdashX*Kbar_Mat*NdashX.transpose()*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                    Mlocal_arr[i] = Mlocal_arr[i] + N.transpose()*N*rhoC*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                    //For Flocal
                    Eigen::MatrixXf Flocal(NumNodesPerElem,1);
                    //Flocal_arr[i] = Flocal_arr[i] + Flocal*F*N*N.transpose()*rhoC*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                    Flocal_arr[i] = Flocal_arr[i] + Flocal;     
                } 
            }
        }
    }

    //Assembly
    for (int i=0;i<NumOfElems;i++)
    {
        if(i==99)
        {
            int debug =1;
        }
        int nodeid1 = ElemArrQuad(i,0);
        int nodeid2 = ElemArrQuad(i,1);
        int nodeid3 = ElemArrQuad(i,2);
        int nodeid4 = ElemArrQuad(i,3);
        Eigen::MatrixXd nodeIds(1, 4);
        nodeIds << nodeid1, nodeid2, nodeid3, nodeid4; //comma initialization
        //std::cout << "nodeIds" << nodeIds << std::endl;
       
        for (int m=0;m<NumNodesPerElem;m++)
        {   
            for (int n=0;n<NumNodesPerElem;n++)
            {
                KGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) = KGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) + Klocal_arr[i](m,n);
                MGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) = MGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) + Mlocal_arr[i](m,n);
            }
            FGlobal(int(nodeIds(0,m)),0) = FGlobal(int(nodeIds(0,m)),0) + Flocal_arr[i](m,0);
        }
    }
   //std::cout << "KGlobal " << KGlobal << std::endl;    
    return 0;
}

int GetTriaKMandF2DFE(int NumOfElems, Eigen::MatrixXf GNodeCordsArr, Eigen::MatrixXd  ElemArrQuad, int order, int bondaryCondition, float F, Eigen::MatrixXf  &KGlobal, Eigen::MatrixXf  &MGlobal, Eigen::MatrixXf  &FGlobal)
{
    
    int NElemsH =20;
    int NElemsV =20;
    int NPerElems =4;
    int dim = 2;
    float xMin = 0.0;
    float yMin = 0.0;
    float xMax = 1.0;
    float yMax = 1.0;
    double kbar = 385;
    double rho = 6000;
    float rhoC = 3.8151e6; 
    Eigen::MatrixXf Kbar_Mat = Eigen::MatrixXf::Identity(dim,dim);
    Kbar_Mat = kbar*Kbar_Mat;

   
    int NumNodesTotal=(NElemsH+1)*(NElemsV+1);
    Eigen::MatrixXf Klocal_arr[800] = {}; //initialize to NumOfElems
    Eigen::MatrixXf Mlocal_arr[800] = {}; //initialize to NumOfElems
    Eigen::MatrixXf Flocal_arr[800] = {};
    int NumNodesPerElem = 3;
    if (order==1) 
        NumNodesPerElem=3;  
    else if (order==2) 
         NumNodesPerElem=6;

    float quadPts[3] = {0.1666, 0.166, 0.666};   //-1/sqrt(3),1/sqrt(3)  
    int NumQuadPts = sizeof(quadPts)/sizeof(int);
    double WeightquadPts[3] = {0.33,0.33,0.33};   
    
    double L = 1.0;
    double Le = L/NumOfElems; 
    float zDashx = 2/Le;

    for (int i=0;i<NumOfElems;i++)
    {
        Klocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,NumNodesPerElem);
        Mlocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,NumNodesPerElem);
        Flocal_arr[i] = Eigen::MatrixXf::Zero(NumNodesPerElem,1);
    }
    

    for (int i=0;i<NumOfElems;i++)
    {
        for (int j=0;j<NumQuadPts;j++)
        {
            for (int k=0;k<NumQuadPts;k++)
            {
                //For Klocal
                Eigen::MatrixXf Klocal(NumNodesPerElem,NumNodesPerElem);
                if (order == 1)
                {
                    Eigen::MatrixXf N = GetShapeFn(quadPts[j],quadPts[k],1);
                    Eigen::MatrixXf NwZ = GetShapeFnGradientwZ(quadPts[j],quadPts[k],1);
                    Eigen::MatrixXf NwEta = GetShapeFnGradientwEta(quadPts[j],quadPts[k],1);

                    int nodeid1 = ElemArrQuad(i,0);
                    int nodeid2 = ElemArrQuad(i,1);
                    int nodeid3 = ElemArrQuad(i,2);

                    Eigen::MatrixXf nodeIds(1, 4);
                    nodeIds << nodeid1, nodeid2, nodeid3; //comma initialization
                    //std::cout << "nodeIds" << nodeIds << std::endl;

                    Eigen::MatrixXf xCoords(3, 1);
                    xCoords << GNodeCordsArr(nodeid1,0), GNodeCordsArr(nodeid2,0), GNodeCordsArr(nodeid3,0); //comma initialization
                    //std::cout << "xCoords" << xCoords << std::endl;
            
                     Eigen::MatrixXf yCoords(3, 1);
                    yCoords << GNodeCordsArr(nodeid1,1), GNodeCordsArr(nodeid2,1), GNodeCordsArr(nodeid3,1); //comma initialization
                    //std::cout << "yCoords" << yCoords << std::endl;

                    Eigen::MatrixXf Jacobian = Eigen::MatrixXf::Zero(dim,dim);
                    Jacobian << NwZ*xCoords, NwEta*xCoords, NwZ*yCoords, NwEta*yCoords;
                    //std::cout << "Jacobian" << Jacobian << std::endl;
                    
                    Eigen::MatrixXf invJacobian = Eigen::MatrixXf::Zero(dim,dim);
                    invJacobian = Jacobian.inverse();
                    //std::cout << "invJacobian" << invJacobian << std::endl;

                    Eigen::MatrixXf NdashZEta = Eigen::MatrixXf::Zero(NumNodesPerElem,dim);
                    NdashZEta.col(0) << NwZ(0,0),NwZ(0,1),NwZ(0,2),NwZ(0,3);
                    NdashZEta.col(1) << NwEta(0,0),NwEta(0,1),NwEta(0,2),NwEta(0,3);
                    //std::cout << "NdashZEta" << NdashZEta << std::endl;

                    Eigen::MatrixXf NdashX = Eigen::MatrixXf::Zero(NumNodesPerElem,dim);
                    NdashX = NdashZEta*invJacobian;
                    // std::cout << "NdashX" << NdashX << std::endl;
                    // std::cout << "Kbar_Mat" << Kbar_Mat << std::endl;

                    Klocal_arr[i] = Klocal_arr[i] + NdashX*Kbar_Mat*NdashX.transpose()*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                    Mlocal_arr[i] = Mlocal_arr[i] + N.transpose()*N*rhoC*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                    //For Flocal
                    Eigen::MatrixXf Flocal(NumNodesPerElem,1);
                    //Flocal_arr[i] = Flocal_arr[i] + Flocal*F*N*N.transpose()*rhoC*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                    Flocal_arr[i] = Flocal_arr[i] + Flocal;     
                } 
            }
        }
    }

    //Assembly
    for (int i=0;i<NumOfElems;i++)
    {
        if(i==99)
        {
            int debug =1;
        }
        int nodeid1 = ElemArrQuad(i,0);
        int nodeid2 = ElemArrQuad(i,1);
        int nodeid3 = ElemArrQuad(i,2);
        Eigen::MatrixXd nodeIds(1, 3);
        nodeIds << nodeid1, nodeid2, nodeid3; //comma initialization
        //std::cout << "nodeIds" << nodeIds << std::endl;
       
        for (int m=0;m<NumNodesPerElem;m++)
        {   
            for (int n=0;n<NumNodesPerElem;n++)
            {
                KGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) = KGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) + Klocal_arr[i](m,n);
                MGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) = MGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) + Mlocal_arr[i](m,n);
            }
            FGlobal(int(nodeIds(0,m)),0) = FGlobal(int(nodeIds(0,m)),0) + Flocal_arr[i](m,0);
        }
    }
   //std::cout << "KGlobal " << KGlobal << std::endl;    
    return 0;
}

Eigen::MatrixXf GetDGlobalQuadElem(int NElemsH, int NElemsV, int order, int dim, int bondaryCondition, float F, int timesteps, int alpha)
{
    int NumNodesTotal=(NElemsH+1)*(NElemsV+1);
    int NumOfElems = NElemsH*NElemsV;
    int NPerElems =4;
    float xMin = 0.0;
    float yMin = 0.0;
    float xMax = 1.0;
    float yMax = 1.0;

    float deltaT = 2.085;

    if (order==1) 
    {
        NumNodesTotal = (NElemsH+1)*(NElemsV+1);//linear
        NPerElems = 4;
    }
    else if (order==2) 
    {
        NumNodesTotal = (2*NElemsH+1)*(2*NElemsV+1);//Quadratic
        NPerElems = 8;
    }
    
    //createmesh
    Eigen::MatrixXd  ElemArrQuad = GetElemArr(NElemsH, NElemsV, NPerElems);
    Eigen::MatrixXf GNodeCordsArr = GetGNodeCordsArr(NElemsH, NElemsV, NPerElems ,dim, xMin, yMin, xMax, yMax);
  
    Eigen::MatrixXf  KGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf  MGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf  FGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf  DGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);

    //Getting KGlobal, MGlobal, FGlobal
   GetQuadKMandF2DFE(NumOfElems, GNodeCordsArr, ElemArrQuad, order, bondaryCondition, F, KGlobal, MGlobal, FGlobal);

    int nDirichletnodes = NElemsH+NElemsH+2;
    int nNotDirchiletnodes = NumNodesTotal - nDirichletnodes;

    //Applying Dirichilet Boundary conditions

//      //Apply Dirichlet's BCs
//     /* for the following sets of boundary conditions and forcing functions:
//        1. (2D Quadrilateral Mesh): (x ∈ [0, 1], y ∈ [0, 1], use a 20 x 20 element mesh.)
//        g = 300 K along x = 0 m (left edge) and g = 310 K along x = 1.0 m (right edge)
//        Initial condition: u(x, y, 0) = 300K for x < 0.5m and u(x, y, 0) = 300 + 20(x − 0.5)K for x ≥ 0.5m

//     where  f_tilde and f_bar are constant values.*/
//     //κij = ¯κδij , where ¯κ = 385 watt.m−1K−1, ρ = 6000 kg.m−3 and ρc = 3.8151 x 106 N.m−2K−1

//     std::ofstream out ("output.txt");
    //std::cout << "I wrote data to the file\n";
    // ifstream is just like cin:

    //g = 300 K along x = 0 m (left edge) and g = 310 K along x = 1.0 m (right edge)
    Eigen::MatrixXf  DirichletBC = Eigen::MatrixXf::Zero(nDirichletnodes,2);
    //left edge
    for (int i=0;i<NElemsH+1;i++)
    {
        DirichletBC(i,0) = i*NElemsH+i;
        DirichletBC(i,1) = 300;
    }
    //right edge
    for (int i=0;i<NElemsH+1;i++)
    {
        DirichletBC(i+NElemsH+1,0) = NElemsH+((NElemsH+1)*i);
        DirichletBC(i+NElemsH+1,1) = 310;
    }

    //KDGlobal, MDGlobal and FDGlobal
    Eigen::MatrixXf  KDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,nNotDirchiletnodes);
    Eigen::MatrixXf  MDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,nNotDirchiletnodes);
    Eigen::MatrixXf  FDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    Eigen::MatrixXf  NotDirchNodes = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    
    //Non Dirchlet Nodes matrix
    int count = 0;
    for (int i=0;i<NElemsH+1;i++)
    {
        for (int j=0;j<NElemsH+1;j++)
        {
            if(j!=0 && j!=NElemsH)
            {
                NotDirchNodes(count,0) = i*(NElemsH+1)+j;
                count++;
            }
        }
        
    }

    //KDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nNotDirchiletnodes;j++)
        {
            KDGlobal(i,j)=KGlobal(int(NotDirchNodes(i,0)),int(NotDirchNodes(j,0)));
        }
    }

    // MDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nNotDirchiletnodes;j++)
        {
            MDGlobal(i,j)=MGlobal(int(NotDirchNodes(i,0)),int(NotDirchNodes(j,0)));
        }
    }

    Eigen::MatrixXf  F_forSubstract = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nDirichletnodes;j++)
        {
            F_forSubstract(i,0)=KGlobal(int(NotDirchNodes(i,0)),int(DirichletBC(j,0)))*DirichletBC(j,1) + F_forSubstract(i,0);
        }
    }

    // FDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        FDGlobal(i,0)=FGlobal(int(NotDirchNodes(i,0)),0) - F_forSubstract(i,0);
    }

    // time discretization
    //  Initial condition: u(x, y, 0) = 300K for x < 0.5m and u(x, y, 0) = 300 + 20(x − 0.5)K for x ≥ 0.5m
    Eigen::MatrixXf  UValues = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    for (int i=0;i<NumNodesTotal;i++)
    {
        if(GNodeCordsArr(i,0)<0.5)
        {
            UValues(i,0)=300;
        }
        if(GNodeCordsArr(i,0)>=0.5)
        {
            UValues(i,0)=300+10*(GNodeCordsArr(i,0)-0.5);
        }
        DGlobal(i,0) = UValues(i,0);
    }

    //Applying Vmethod
    Eigen::MatrixXf  Dn = Eigen::MatrixXf::Zero(nNotDirchiletnodes,timesteps+1);
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        int nodeid = NotDirchNodes(i,0);
        Dn(i,0) = DGlobal(nodeid,0);
    }

    Eigen::MatrixXf  Vn = Eigen::MatrixXf::Zero(nNotDirchiletnodes,timesteps+1);
    Vn.col(0) << MDGlobal.inverse()*(FDGlobal - KDGlobal*Dn.col(0));

    Eigen::MatrixXf  dTilde = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    for (int i=0;i<timesteps;i++)
    { 
        dTilde = Dn.col(i) + (1-alpha)*deltaT*Vn.col(i);
        Vn.col(i+1) = (MDGlobal + alpha*deltaT*KDGlobal).inverse()*(FDGlobal - KDGlobal*dTilde);
        Dn.col(i+1) = dTilde + alpha*deltaT*Vn.col(i+1);
    }

    for (int i=0;i<nDirichletnodes;i++)
    { 
        DGlobal(int(DirichletBC(i,0)),0)=int(DirichletBC(i,1));
    }
    for (int i=0;i<nNotDirchiletnodes;i++)
    { 
        DGlobal(int(NotDirchNodes(i,0)),0)=int(Dn(i,timesteps-1));
    }
    //std::cout << "DGlobal" << DGlobal << "\n";

    return DGlobal;
}


Eigen::MatrixXf GetDGlobalTriaElem(int NElemsH, int NElemsV, int order, int dim, int bondaryCondition, float F, int timesteps, int alpha)
{
    int NumNodesTotal=(NElemsH+1)*(NElemsV+1);
    int NumOfElems = NElemsH*NElemsV;
    int NPerElems =3;
    float xMin = 0.0;
    float yMin = 0.0;
    float xMax = 1.0;
    float yMax = 1.0;

    float deltaT = 2.085;

    if (order==1) 
    {
        NumNodesTotal = (NElemsH+1)*(NElemsV+1);//linear
        NPerElems = 4;
    }
    else if (order==2) 
    {
        NumNodesTotal = (2*NElemsH+1)*(2*NElemsV+1);//Quadratic
        NPerElems = 8;
    }
    
    //createmesh
    Eigen::MatrixXd  ElemArrQuad = GetTriaElemArr(NElemsH, NElemsV, NPerElems);
    Eigen::MatrixXf GNodeCordsArr = GetTriaGNodeCordsArr(NElemsH, NElemsV, NPerElems ,dim, xMin, yMin, xMax, yMax);
  
    Eigen::MatrixXf  KGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf  MGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,NumNodesTotal);
    Eigen::MatrixXf  FGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    Eigen::MatrixXf  DGlobal = Eigen::MatrixXf::Zero(NumNodesTotal,1);

    //Getting KGlobal, MGlobal, FGlobal
   GetTriaKMandF2DFE(NumOfElems, GNodeCordsArr, ElemArrQuad, order, bondaryCondition, F, KGlobal, MGlobal, FGlobal);

    int nDirichletnodes = NElemsH+NElemsH+2;
    int nNotDirchiletnodes = NumNodesTotal - nDirichletnodes;

    //Applying Dirichilet Boundary conditions

//      //Apply Dirichlet's BCs
//     /* for the following sets of boundary conditions and forcing functions:
//        1. (2D Quadrilateral Mesh): (x ∈ [0, 1], y ∈ [0, 1], use a 20 x 20 element mesh.)
//        g = 300 K along x = 0 m (left edge) and g = 310 K along x = 1.0 m (right edge)
//        Initial condition: u(x, y, 0) = 300K for x < 0.5m and u(x, y, 0) = 300 + 20(x − 0.5)K for x ≥ 0.5m

//     where  f_tilde and f_bar are constant values.*/
//     //κij = ¯κδij , where ¯κ = 385 watt.m−1K−1, ρ = 6000 kg.m−3 and ρc = 3.8151 x 106 N.m−2K−1

//     std::ofstream out ("output.txt");
    //std::cout << "I wrote data to the file\n";
    // ifstream is just like cin:

    //g = 300 K along x = 0 m (left edge) and g = 310 K along x = 1.0 m (right edge)
    Eigen::MatrixXf  DirichletBC = Eigen::MatrixXf::Zero(nDirichletnodes,2);
    //left edge
    for (int i=0;i<NElemsH+1;i++)
    {
        DirichletBC(i,0) = i*NElemsH+i;
        DirichletBC(i,1) = 300;
    }
    //right edge
    for (int i=0;i<NElemsH+1;i++)
    {
        DirichletBC(i+NElemsH+1,0) = NElemsH+((NElemsH+1)*i);
        DirichletBC(i+NElemsH+1,1) = 305;
    }

    //KDGlobal, MDGlobal and FDGlobal
    Eigen::MatrixXf  KDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,nNotDirchiletnodes);
    Eigen::MatrixXf  MDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,nNotDirchiletnodes);
    Eigen::MatrixXf  FDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    Eigen::MatrixXf  NotDirchNodes = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    
    //Non Dirchlet Nodes matrix
    int count = 0;
    for (int i=0;i<NElemsH+1;i++)
    {
        for (int j=0;j<NElemsH+1;j++)
        {
            if(j!=0 && j!=NElemsH)
            {
                NotDirchNodes(count,0) = i*(NElemsH+1)+j;
                count++;
            }
        }
        
    }

    //KDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nNotDirchiletnodes;j++)
        {
            KDGlobal(i,j)=KGlobal(int(NotDirchNodes(i,0)),int(NotDirchNodes(j,0)));
        }
    }

    // MDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nNotDirchiletnodes;j++)
        {
            MDGlobal(i,j)=MGlobal(int(NotDirchNodes(i,0)),int(NotDirchNodes(j,0)));
        }
    }

    Eigen::MatrixXf  F_forSubstract = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nDirichletnodes;j++)
        {
            F_forSubstract(i,0)=KGlobal(int(NotDirchNodes(i,0)),int(DirichletBC(j,0)))*DirichletBC(j,1) + F_forSubstract(i,0);
        }
    }

    // FDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        FDGlobal(i,0)=FGlobal(int(NotDirchNodes(i,0)),0) - F_forSubstract(i,0);
    }

    // time discretization
    //  Initial condition: u(x, y, 0) = 300K for x < 0.5m and u(x, y, 0) = 300 + 20(x − 0.5)K for x ≥ 0.5m
    Eigen::MatrixXf  UValues = Eigen::MatrixXf::Zero(NumNodesTotal,1);
    for (int i=0;i<NumNodesTotal;i++)
    {
        if(GNodeCordsArr(i,0)<0.5)
        {
            UValues(i,0)=300;
        }
        if(GNodeCordsArr(i,0)>=0.5)
        {
            UValues(i,0)=300+10*(GNodeCordsArr(i,0)-0.5);
        }
        DGlobal(i,0) = UValues(i,0);
    }

    //Applying Vmethod
    Eigen::MatrixXf  Dn = Eigen::MatrixXf::Zero(nNotDirchiletnodes,timesteps+1);
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        int nodeid = NotDirchNodes(i,0);
        Dn(i,0) = DGlobal(nodeid,0);
    }

    Eigen::MatrixXf  Vn = Eigen::MatrixXf::Zero(nNotDirchiletnodes,timesteps+1);
    Vn.col(0) << MDGlobal.inverse()*(FDGlobal - KDGlobal*Dn.col(0));

    Eigen::MatrixXf  dTilde = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    for (int i=0;i<timesteps;i++)
    { 
        dTilde = Dn.col(i) + (1-alpha)*deltaT*Vn.col(i);
        Vn.col(i+1) = (MDGlobal + alpha*deltaT*KDGlobal).inverse()*(FDGlobal - KDGlobal*dTilde);
        Dn.col(i+1) = dTilde + alpha*deltaT*Vn.col(i+1);
    }

    for (int i=0;i<nDirichletnodes;i++)
    { 
        DGlobal(int(DirichletBC(i,0)),0)=int(DirichletBC(i,1));
    }
    for (int i=0;i<nNotDirchiletnodes;i++)
    { 
        DGlobal(int(NotDirchNodes(i,0)),0)=int(Dn(i,timesteps-1));
    }
    //std::cout << "DGlobal" << DGlobal << "\n";

    return DGlobal;
}

int main()
{
    {
    int NElemsH =20;//horizontal
    int NElemsV =20;// vertical
    int NumOfElems = NElemsH*NElemsV;
    int NPerElems =4;
    int dim = 2;
    float xMin = 0.0;
    float yMin = 0.0;
    float xMax = 1.0;
    float yMax = 1.0;

    int order = 1;
    //Apply Dirichlet's BCs
    //     /* for the following sets of boundary conditions and forcing functions:
    //        1. (2D Quadrilateral Mesh): (x ∈ [0, 1], y ∈ [0, 1], use a 20 x 20 element mesh.)
    //        g = 300 K along x = 0 m (left edge) and g = 310 K along x = 1.0 m (right edge)
    //        Initial condition: u(x, y, 0) = 300K for x < 0.5m and u(x, y, 0) = 300 + 20(x − 0.5)K for x ≥ 0.5m

    //     where  f_tilde and f_bar are constant values.*/
    //     //κij = ¯κδij , where ¯κ = 385 watt.m−1K−1, ρ = 6000 kg.m−3 and ρc = 3.8151 x 106 N.m−2K−1

    std::ofstream out ("ElemArrHW3.txt");
    std::cout << "write data to the file\n";
    //createmesh for output
    Eigen::MatrixXd  ElemArrQuad = GetElemArr(NElemsH, NElemsV, NPerElems);
    std::ofstream out1 ("NodeCoordsHW3.txt");
     Eigen::MatrixXf GNodeCordsArr = GetGNodeCordsArr(NElemsH, NElemsV, NPerElems ,dim, xMin, yMin, xMax, yMax);

    for(int i=0; i< NElemsH*NElemsV;i++)
    {
        int k=0;
        int nGNodeId1 = ElemArrQuad(i,k) ;
        int nGNodeId2 = ElemArrQuad(i,k+1);
        int nGNodeId3 = ElemArrQuad(i,k+2);
        int nGNodeId4 = ElemArrQuad(i,k+3);

        ElemArrQuad(i,k)  = nGNodeId1+1;
        ElemArrQuad(i,k+1) = nGNodeId2+1;
        ElemArrQuad(i,k+2) = nGNodeId3+1;
        ElemArrQuad(i,k+3) = nGNodeId4+1;
    }
    std::cout << "ElemArrQuad" << ElemArrQuad << std::endl;

    out << "ElemArrQuad" "\n";
    out << ElemArrQuad << "\n";
    out.close();

    out1 << "GNodeCordsArr" "\n";
    out1 << GNodeCordsArr << "\n";
    out1.close();

    int bondaryCondition =1;
    float F = 0;
    int timesteps = 10;
    double alpha = 1;
    Eigen::MatrixXf  DGlobal = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal << std::endl;
    std::ofstream out2 ("DGlobalHW3BC1.txt");
    out2 << "DGlobal" "\n";
    out2 << DGlobal << "\n";
    out2.close();
    
//     std::cout << "Boundary Condition 2" << std::endl;
    timesteps = 25;
    alpha = 1;
    Eigen::MatrixXf  DGlobal2 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    std::ofstream out3 ("DGlobalHW3BC2.txt");
    out3 << "DGlobal" "\n";
    out3 << DGlobal2 << "\n";
    out3.close();
    //std::cout << "DGlobal" << DGlobal2 << std::endl;
//     out << "DGlobalBC2" << DGlobalBC2 << "\n";

//     std::cout << "Boundary Condition 3" << std::endl;
    timesteps = 100;
    alpha = 1;
    Eigen::MatrixXf  DGlobal3 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal3 << std::endl;
    std::ofstream out4 ("DGlobalHW3BC3.txt");
    out4 << "DGlobal" "\n";
    out4 << DGlobal3 << "\n";
    out4.close();
//out.close();

  timesteps = 10;
     alpha = 0.5;
    Eigen::MatrixXf  DGlobal4 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal4 << std::endl;
   
     std::ofstream out5 ("DGlobalHW3BC4.txt");
    out5 << "DGlobal" "\n";
    out5 << DGlobal4 << "\n";
    out5.close();
    
//     std::cout << "Boundary Condition 2" << std::endl;
    timesteps = 25;
    alpha = 0.5;
    Eigen::MatrixXf  DGlobal5 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal5 << std::endl;
//     out << "DGlobalBC2" << DGlobalBC2 << "\n";
    std::ofstream out6 ("DGlobalHW3BC5.txt");
    out6 << "DGlobal" "\n";
    out6 << DGlobal5 << "\n";
    out6.close();

//     std::cout << "Boundary Condition 3" << std::endl;
    timesteps = 100;
    alpha = 0;
    Eigen::MatrixXf  DGlobal6 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal6 << std::endl;
std::ofstream out7 ("DGlobalHW3BC6.txt");
    out7 << "DGlobal" "\n";
    out7 << DGlobal6 << "\n";
    out7.close();


  timesteps = 10;
     alpha = 0;
    Eigen::MatrixXf  DGlobal7 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal7 << std::endl;
    out << "DGlobal" "\n";
    out << DGlobal << "\n";
    std::ofstream out8 ("DGlobalHW3BC7.txt");
    out8 << "DGlobal" "\n";
    out8 << DGlobal7 << "\n";
    out8.close();

//     std::cout << "Boundary Condition 2" << std::endl;
    timesteps = 25;
    alpha = 0;
    Eigen::MatrixXf  DGlobal9 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal9 << std::endl;
//     out << "DGlobalBC2" << DGlobalBC2 << "\n";
std::ofstream out9 ("DGlobalHW3BC8.txt");
    out9 << "DGlobal" "\n";
    out9 << DGlobal9 << "\n";
    out9.close();

//     std::cout << "Boundary Condition 3" << std::endl;
    timesteps = 100;
    alpha = 0.5;
    Eigen::MatrixXf  DGlobal8 = GetDGlobalQuadElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal8 << std::endl;
std::ofstream out10 ("DGlobalHW3BC9.txt");
    out10 << "DGlobal" "\n";
    out10 << DGlobal8 << "\n";
    out10.close();

//out.close();
   
    }
    {
    //For triangular elems
    int NElemsH =40;//horizontal
    int NElemsV =40;// vertical
    int NumOfElems = NElemsH*NElemsV;
    int NPerElems =3;
    int dim = 2;
    float xMin = 0.0;
    float yMin = 0.0;
    float xMax = 1.0;
    float yMax = 1.0;
int bondaryCondition =1;
    float F = 0;

    int order = 1;
    Eigen::MatrixXd  ElemArrTria = GetTriaElemArr(NElemsH, NElemsV, NPerElems);
     //std::cout << "ElemArrTria" << ElemArrTria << std::endl;
    Eigen::MatrixXf GNodeCordsArr = GetTriaGNodeCordsArr(NElemsH, NElemsV, NPerElems ,dim, xMin, yMin, xMax, yMax);
     //std::cout << "GetTriaGNodeCordsArr" << GNodeCordsArr << std::endl;

     int timesteps = 10;
    int alpha = 1;
    Eigen::MatrixXf  DGlobal = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal << std::endl;
    //out << "DGlobal" "\n";
    //out << DGlobal << "\n";
    
//     std::cout << "Boundary Condition 2" << std::endl;
    timesteps = 25;
    alpha = 1;
    Eigen::MatrixXf  DGlobal2 = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal2 << std::endl;
//     out << "DGlobalBC2" << DGlobalBC2 << "\n";

//     std::cout << "Boundary Condition 3" << std::endl;
    timesteps = 100;
    alpha = 1;
    Eigen::MatrixXf  DGlobal3 = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    std::cout << "DGlobal3" << DGlobal3 << std::endl;
//out.close();

int timesteps = 10;
    int alpha = 0.5;
    Eigen::MatrixXf  DGlobal = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal << std::endl;
    //out << "DGlobal" "\n";
    //out << DGlobal << "\n";
    
//     std::cout << "Boundary Condition 2" << std::endl;
    timesteps = 25;
    alpha = 0.5;
    Eigen::MatrixXf  DGlobal2 = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal2 << std::endl;
//     out << "DGlobalBC2" << DGlobalBC2 << "\n";

//     std::cout << "Boundary Condition 3" << std::endl;
    timesteps = 100;
    alpha = 0.5;
    Eigen::MatrixXf  DGlobal3 = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    std::cout << "DGlobal3" << DGlobal3 << std::endl;

    int timesteps = 10;
    int alpha = 0.5;
    Eigen::MatrixXf  DGlobal = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal << std::endl;
    //out << "DGlobal" "\n";
    //out << DGlobal << "\n";
    
//     std::cout << "Boundary Condition 2" << std::endl;
    timesteps = 25;
    alpha = 0.5;
    Eigen::MatrixXf  DGlobal2 = GetDGlobalTriaElem(NElemsH, NElemsV, order, dim, bondaryCondition, F,timesteps, alpha);
    //std::cout << "DGlobal" << DGlobal2 << std::endl;
//     out << "DGlobalBC2" << DGlobalBC2 << "\n";   

//     Eigen::MatrixXf GNodeCordsArr = GetTriaGNodeCordsArr(NElemsH, NElemsV, NPerElems ,dim, xMin, yMin, xMax, yMax);
     }
    
    return 0;
}
