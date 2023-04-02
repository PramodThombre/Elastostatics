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

//Linear Shape functions and its derivatives

Eigen::MatrixXf GetShapeFn(double z, double eta, double k, int order)
{
    Eigen::MatrixXf N_ShapeFn = Eigen::MatrixXf::Zero(1,8);
    N_ShapeFn(0,0) = (1-z)*(1-eta)*(1-k)/8;
    N_ShapeFn(0,1) = (1+z)*(1-eta)*(1-k)/8;
    N_ShapeFn(0,2) = (1+z)*(1+eta)*(1-k)/8;
    N_ShapeFn(0,3) = (1-z)*(1+eta)*(1-k)/8;

    N_ShapeFn(0,4) = (1-z)*(1-eta)*(1+k)/8;
    N_ShapeFn(0,5) = (1+z)*(1-eta)*(1+k)/8;
    N_ShapeFn(0,6) = (1+z)*(1+eta)*(1+k)/8;
    N_ShapeFn(0,7) = (1-z)*(1+eta)*(1+k)/8;
    return N_ShapeFn;
}
Eigen::MatrixXf GetShapeFnGradientwZ(double z, double eta, double k, int order)
{
    Eigen::MatrixXf N_ShapeFnGradient = Eigen::MatrixXf::Zero(1,8);
    N_ShapeFnGradient(0,0) = -(1-eta)*(1-k)/8;
    N_ShapeFnGradient(0,1) = (1-eta)*(1-k)/8;
    N_ShapeFnGradient(0,2) = (1+eta)*(1-k)/8;
    N_ShapeFnGradient(0,3) = -(1+eta)*(1-k)/8;
    N_ShapeFnGradient(0,4) = -(1-eta)*(1+k)/8;
    N_ShapeFnGradient(0,5) = (1-eta)*(1+k)/8;
    N_ShapeFnGradient(0,6) = (1+eta)*(1+k)/8;
    N_ShapeFnGradient(0,7) = -(1+eta)*(1+k)/8;

    return N_ShapeFnGradient;
}
Eigen::MatrixXf GetShapeFnGradientwEta(double z, double eta, double k, int order)
{
    Eigen::MatrixXf N_ShapeFnGradient = Eigen::MatrixXf::Zero(1,8);
    N_ShapeFnGradient(0,0) = -(1-z)*(1-k)/8;
    N_ShapeFnGradient(0,1) = -(1+z)*(1-k)/8;
    N_ShapeFnGradient(0,2) = (1+z)*(1-k)/8;
    N_ShapeFnGradient(0,3) = (1-z)*(1-k)/8;
    N_ShapeFnGradient(0,4) = -(1-z)*(1+k)/8;
    N_ShapeFnGradient(0,5) = -(1+z)*(1+k)/8;
    N_ShapeFnGradient(0,6) = (1+z)*(1+k)/8;
    N_ShapeFnGradient(0,7) = (1-z)*(1+k)/8;
    return N_ShapeFnGradient;
}
Eigen::MatrixXf GetShapeFnGradientwK(double z, double eta, double k, int order)
{
    Eigen::MatrixXf N_ShapeFnGradient = Eigen::MatrixXf::Zero(1,8);
    N_ShapeFnGradient(0,0) = -(1-z)*(1-eta)/8;;
    N_ShapeFnGradient(0,1) = -(1+z)*(1-eta)/8;;
    N_ShapeFnGradient(0,2) = -(1+z)*(1+eta)/8;;
    N_ShapeFnGradient(0,3) = -(1-z)*(1+eta)/8;;
    N_ShapeFnGradient(0,4) = (1-z)*(1-eta)/8;
    N_ShapeFnGradient(0,5) = (1+z)*(1-eta)/8;
    N_ShapeFnGradient(0,6) = (1+z)*(1+eta)/8;
    N_ShapeFnGradient(0,7) = (1-z)*(1+eta)/8;
    return N_ShapeFnGradient;
}

Eigen::MatrixXd GetElemArr(int NXElems, int NYElems, int NZElems, int NNodesPerElem)
{
    //Elem_LocaltoGlobalNodeArr(ElemId,LocalNodeId) = GlobalNodeId;
    Eigen::MatrixXd Elem_LocaltoGlobalNodeArr = Eigen::MatrixXd::Zero(NXElems*NYElems*NZElems,NNodesPerElem);
    int count = 0;
    for(int i=0;i<NZElems;i++)//in z direction
    {
        //std::cout << "k," << k << std::endl;
        for(int j=0;j<NYElems;j++)//in y direction
        {
           //std::cout << "i," << i << std::endl;
            for(int k=0;k<NXElems;k++)//in x direction
            {
                //first looping through all Z dir elems then to Y dir and then to X dir
                //std::cout << "k:" << k << ",j:" << j <<",i:" <<i << std::endl;

                int nElemId = k+(j)*NXElems+(i)*NXElems*NYElems;
                int nGNodeId1 =nElemId+(j)+(i)*(NXElems+NYElems+1);
                int nGNodeId2 =nElemId+(j+1)+(i)*(NXElems+NYElems+1);
                int nGNodeId3 =nElemId+(j+1)+NXElems+1+(i)*(NXElems+NYElems+1);
                int nGNodeId4 =nElemId+(j+1)+NXElems+(i)*(NXElems+NYElems+1);
                
                int nGNodeId5 =nElemId+(j)+(i+1)*(NXElems+NYElems+1)+NXElems*NYElems;
                int nGNodeId6 =nElemId+(j+1)+(i+1)*(NXElems+NYElems+1)+NXElems*NYElems;
                int nGNodeId7 =nElemId+(j+1)+NXElems+1+(i+1)*(NXElems+NYElems+1)+NXElems*NYElems;
                int nGNodeId8 =nElemId+(j+1)+NXElems+(i+1)*(NXElems+NYElems+1)+NXElems*NYElems;

                int nPos=0;
                std::cout << "nElemId," << nElemId << std::endl;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos) = nGNodeId1;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+1) = nGNodeId2;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+2) = nGNodeId3;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+3) = nGNodeId4;
                std::cout << "nElemId," << nElemId << " nGNodeId1," << nGNodeId1 << std::endl;
                std::cout << "nElemId," << nElemId << " nGNodeId2," << nGNodeId2 << std::endl;
                std::cout << "nElemId," << nElemId << " nGNodeId3," << nGNodeId3<< std::endl;
                std::cout << "nElemId," << nElemId << " nGNodeId4," << nGNodeId4 << std::endl;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+4) = nGNodeId5;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+5) = nGNodeId6;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+6) = nGNodeId7;
                Elem_LocaltoGlobalNodeArr(nElemId,nPos+7) = nGNodeId8;
                std::cout << "nElemId," << nElemId << " nGNodeId5," << nGNodeId5 << std::endl;
                std::cout << "nElemId," << nElemId << " nGNodeId6," << nGNodeId6 << std::endl;
                std::cout << "nElemId," << nElemId << " nGNodeId7," << nGNodeId7<< std::endl;
                std::cout << "nElemId," << nElemId << " nGNodeId8," << nGNodeId8 << std::endl;
            
                count++;
                //std::cout << "count," << count << std::endl;
            }
        }
    }
    return Elem_LocaltoGlobalNodeArr;
}

Eigen::MatrixXf GetGNodeCordsArr(int NXElems, int NYElems, int NZElems,  int NNodesPerElem, int dim, float xMin, float yMin, float zMin,float xMax, float yMax, float zMax)
{
    float xLength = xMax - xMin;
    float yLength = yMax - yMin;
    float zLength = zMax - zMin;
    float xLe = xLength/NXElems; //length of each elem in x direction
    float yLe = yLength/NYElems; //length of each elem in y direction
    float zLe = zLength/NZElems; //length of each elem in z direction
    int nNodesX = NXElems+1;
    int nNodesY = NYElems+1;
    int nNodesZ = NZElems+1;

    Eigen::MatrixXf Elem_GlobalNodeCordsArr = Eigen::MatrixXf::Zero((NXElems+1)*(NYElems+1)*(NZElems+1),dim);
    int count = 0;

    for(int i=0;i<NZElems+1;i++)//in z direction
    {
        //std::cout << "k," << k << std::endl;
        for(int j=0;j<NYElems+1;j++)//in y direction
        {
           //std::cout << "i," << i << std::endl;
            for(int k=0;k<NXElems+1;k++)//in x direction
            {
            //std::cout << "k:" << k << ",j:" << j <<",i:" <<i << std::endl;
            int nodeID = k+(j)*(NXElems+1)+(i)*(NXElems+1)*(NYElems+1);
            std::cout << "nodeID," << nodeID << std::endl;

                Elem_GlobalNodeCordsArr(nodeID,0)=(k)*xLe; 
                Elem_GlobalNodeCordsArr(nodeID,1)=(j)*yLe; 
                Elem_GlobalNodeCordsArr(nodeID,2)=(i)*zLe;   

                 count++;
                //std::cout << "count," << count << std::endl;
            }
        }
    }
    std::cout << "Elem_GlobalNodeCordsArr," << Elem_GlobalNodeCordsArr << std::endl;
    return Elem_GlobalNodeCordsArr;
}

// Elasticity tensor
double Cijkl(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    double E = 1000 ;//Pa
    double Nu = 0.3 ;
    float lambda = E/(2*(1+Nu));
    float Mu = E*Nu/((1+Nu)*(1-2*Nu));

    return lambda * (i==j) * (k==l) + 2 * Mu * ((i==k)*(j==l) + (i==l)* (j==k));
}

Eigen::MatrixXf GetNodalStiffness(Eigen::MatrixXf NdashX, int dim1, int dim2)
{
    int dim=3;

    Eigen::MatrixXf KAB = Eigen::MatrixXf::Zero(dim,dim);
    Eigen::MatrixXf NdashX_Transpose=NdashX.transpose();

    for (int i=0;i<dim;i++)
    {
        for (int j=0;j<dim;j++)
        {
            for (int k=0;k<dim;k++)
            {
                for (int l=0;l<dim;l++)
                {
                    // std::cout << "KAB(i,k)," << KAB(i,k) << std::endl;
                    // std::cout << "Cijkl(i,j,k,l)" << Cijkl(i,j,k,l) << std::endl;
                    // std::cout << "NdashX_Transpose(l,dim2)" << NdashX_Transpose(l,dim2) << std::endl;
                    KAB(i,k)= KAB(i,k)+ NdashX(dim1,j)*Cijkl(i,j,k,l)*NdashX_Transpose(l,dim2);
                }
            }
        }
    }

    return KAB;
}



int GetKMandF3DFE(int NElemsX, int NElemsY, int NElemsZ, Eigen::MatrixXf GNodeCordsArr, Eigen::MatrixXd  ElemArrQuad, int order, int bondaryCondition, float F, Eigen::MatrixXf  &KGlobal, Eigen::MatrixXf  &MGlobal, Eigen::MatrixXf  &FGlobal)
{
    
    //Inputs initialization
    double E = 1000 ;//Pa
    double Nu = 0.3 ;
    double g1 = 0;//, 
    double g2 = 0.05;// m, 

    int NumNodesTotal = (NElemsX+1)*(NElemsY+1)*(NElemsZ+1);

    float lambda = E/(2*(1+Nu));
    float Mu = E*Nu/((1+Nu)*(1-2*Nu));

    int NumOfElems = NElemsX*NElemsY*NElemsZ;
    int NPerElems =8; //nodes per Elem
    int dim = 3;

    float xMin = 0.0;
    float yMin = 0.0;
    float zMin = 0.0;
    float xMax = 0.1;
    float yMax = 0.1;
    float ZMax = 1.0;
    int numElems = NElemsX*NElemsY*NElemsZ;
    // Eigen::MatrixXf Klocal_arr[numElems] = {}; //initialize to NumOfElems NElemsX*NElemsY*NElemsZ
    // Eigen::MatrixXf Mlocal_arr[numElems] = {}; //initialize to NumOfElems
    // Eigen::MatrixXf Flocal_arr[numElems] = {};
    int NumNodesPerElem = 8;
    if (order==1) 
        NumNodesPerElem=8;  
    else if (order==2) 
        NumNodesPerElem=20;

    double quadPts[3] = {-0.7746,0,0.7746};   //-sqrt(3/5),0,sqrt(3/5)  
    int NumQuadPts = sizeof(quadPts)/sizeof(double);
    double WeightquadPts[3] = {0.5555,0.8888,0.5555};   //5/9,8/9,5/9  

    for (int i=0;i<NumOfElems;i++)
    {
        Eigen::MatrixXf Klocal = Eigen::MatrixXf::Zero(NumNodesPerElem*dim,NumNodesPerElem*dim);
        Eigen::MatrixXf Mlocal = Eigen::MatrixXf::Zero(NumNodesPerElem*dim,NumNodesPerElem*dim);
        Eigen::MatrixXf Flocal = Eigen::MatrixXf::Zero(NumNodesPerElem*dim,1);
        for (int j=0;j<NumQuadPts;j++)
        {
            for (int k=0;k<NumQuadPts;k++)
            {
                for (int l=0;l<NumQuadPts;l++)
                {
                    //For Klocal
                    //Eigen::MatrixXf Klocal(NumNodesPerElem,NumNodesPerElem);
                    if (order == 1)
                    {
                        Eigen::MatrixXf N = GetShapeFn(quadPts[j],quadPts[k],quadPts[l],1);
                        Eigen::MatrixXf NwZ = GetShapeFnGradientwZ(quadPts[j],quadPts[k],quadPts[l],1);
                        Eigen::MatrixXf NwEta = GetShapeFnGradientwEta(quadPts[j],quadPts[k],quadPts[l],1);
                        Eigen::MatrixXf NwK = GetShapeFnGradientwK(quadPts[j],quadPts[k],quadPts[l],1);

                        int nodeid1 = ElemArrQuad(i,0);
                        int nodeid2 = ElemArrQuad(i,1);
                        int nodeid3 = ElemArrQuad(i,2);
                        int nodeid4 = ElemArrQuad(i,3);
                        int nodeid5 = ElemArrQuad(i,4);
                        int nodeid6 = ElemArrQuad(i,5);
                        int nodeid7 = ElemArrQuad(i,6);
                        int nodeid8 = ElemArrQuad(i,7);

                        Eigen::MatrixXf nodeIds(1, 8);
                        nodeIds << nodeid1, nodeid2, nodeid3, nodeid4, nodeid5, nodeid6, nodeid7, nodeid8; //comma initialization
                        //std::cout << "nodeIds" << nodeIds << std::endl;

                        Eigen::MatrixXf xCoords(8, 1);
                        xCoords << GNodeCordsArr(nodeid1,0), GNodeCordsArr(nodeid2,0), GNodeCordsArr(nodeid3,0), GNodeCordsArr(nodeid4,0),
                         GNodeCordsArr(nodeid5,0), GNodeCordsArr(nodeid6,0), GNodeCordsArr(nodeid7,0), GNodeCordsArr(nodeid8,0); //comma initialization
                        //std::cout << "xCoords" << xCoords << std::endl;
                
                        Eigen::MatrixXf yCoords(8, 1);
                        yCoords << GNodeCordsArr(nodeid1,1), GNodeCordsArr(nodeid2,1), GNodeCordsArr(nodeid3,1), GNodeCordsArr(nodeid4,1), //comma initialization
                        GNodeCordsArr(nodeid5,1), GNodeCordsArr(nodeid6,1), GNodeCordsArr(nodeid7,1), GNodeCordsArr(nodeid8,1); //comma initialization
                        //std::cout << "yCoords" << yCoords << std::endl;
                        
                        Eigen::MatrixXf zCoords(8, 1);
                        zCoords << GNodeCordsArr(nodeid1,2), GNodeCordsArr(nodeid2,2), GNodeCordsArr(nodeid3,2), GNodeCordsArr(nodeid4,2), //comma initialization
                        GNodeCordsArr(nodeid5,2), GNodeCordsArr(nodeid6,2), GNodeCordsArr(nodeid7,2), GNodeCordsArr(nodeid8,2); //comma initialization
                        //std::cout << "yCoords" << yCoords << std::endl;
       
                        Eigen::MatrixXf Jacobian = Eigen::MatrixXf::Zero(dim,dim);
                        Jacobian << NwZ*xCoords,   NwEta*xCoords,  NwK*xCoords, NwZ*yCoords, NwEta*yCoords, NwK*yCoords, NwZ*zCoords, NwEta*zCoords, NwK*zCoords;
                        //std::cout << "Jacobian" << Jacobian << std::endl;
                        
                        Eigen::MatrixXf invJacobian = Eigen::MatrixXf::Zero(dim,dim);
                        invJacobian = Jacobian.inverse();
                        //std::cout << "invJacobian" << invJacobian << std::endl;

                        Eigen::MatrixXf NdashZEta = Eigen::MatrixXf::Zero(NumNodesPerElem,dim);
                        NdashZEta.col(0) << NwZ(0,0),NwZ(0,1),NwZ(0,2),NwZ(0,3),NwZ(0,4),NwZ(0,5),NwZ(0,6),NwZ(0,7);
                        NdashZEta.col(1) << NwEta(0,0),NwEta(0,1),NwEta(0,2),NwEta(0,3),NwEta(0,4),NwEta(0,5),NwEta(0,6),NwEta(0,7);
                        NdashZEta.col(2) << NwK(0,0),NwK(0,1),NwK(0,2),NwK(0,3),NwK(0,4),NwK(0,5),NwK(0,6),NwK(0,7);
                        // std::cout << "N" << N << std::endl;
                        // std::cout << "NwZ" << NwZ << std::endl;
                        // std::cout << "NwEta" << NwEta << std::endl;
                        // std::cout << "NwK" << NwK << std::endl;
                        // std::cout << "NdashZEta" << NdashZEta << std::endl;

                        Eigen::MatrixXf NdashX = Eigen::MatrixXf::Zero(NumNodesPerElem,dim);
                        NdashX = NdashZEta*invJacobian;
                        //std::cout << "NdashX" << NdashX << std::endl;
 
                        Eigen::MatrixXf kL = Eigen::MatrixXf::Zero(NumNodesPerElem*dim,NumNodesPerElem*dim);
                        for (int m=0;m<NumNodesPerElem;m++)
                        {
                            for (int n=0;n<NumNodesPerElem;n++)
                            {
                                Eigen::MatrixXf KAB = Eigen::MatrixXf::Zero(dim,dim);
                                KAB = GetNodalStiffness(NdashX, m, n);
                                //std::cout << "KAB" << KAB << std::endl;

                                kL(3*(m),3*(n))=KAB(0,0);
                                kL(3*(m),3*(n)+1)=KAB(0,1);
                                kL(3*(m),3*(n)+2)=KAB(0,2);
                                kL(3*(m)+1,3*(n))=KAB(1,0);
                                kL(3*(m)+1,3*(n)+1)=KAB(1,1);
                                kL(3*(m)+1,3*(n)+2)=KAB(1,2);
                                kL(3*(m)+2,3*(n))=KAB(2,0);
                                kL(3*(m)+2,3*(n)+1)=KAB(2,1);
                                kL(3*(m)+2,3*(n)+2)=KAB(2,2);

                            }
                        }
                        // std::cout << "kL" << kL << std::endl;
                        // std::cout << "Klocal" << Klocal << std::endl;
                        // std::cout << "WeightquadPts[j]*WeightquadPts[k]*WeightquadPts[l]" << WeightquadPts[j]*WeightquadPts[k]*WeightquadPts[l] << std::endl;
                        // std::cout << "Jacobian.determinant()" << Jacobian.determinant() << std::endl;
                        Eigen::MatrixXf Klocalinloop = kL*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k]*WeightquadPts[l];
                        // std::cout << "Klocalinloop" << Klocalinloop << std::endl;
                        // std::cout << "The vector Klocalinloop is of size " << Klocalinloop.size() << std::endl;
                        // std::cout << "The vector Klocal is of size " << Klocal.size() << std::endl;

                        Klocal = Klocal + Klocalinloop;
                           
                        //For Flocal
                        Eigen::MatrixXf Flocalinloop(NumNodesPerElem*dim,1);
                        //Flocal_arr[i] = Flocal_arr[i] + Flocal*F*N*N.transpose()*rhoC*Jacobian.determinant()*WeightquadPts[j]*WeightquadPts[k];
                        Flocal = Flocal + Flocalinloop;     
                    } 
                }
            }
        }
        //std::cout << "Klocal" << Klocal << std::endl;
        if(i==99)
        {
            int debug =1;
        }
        int nodeid1 = ElemArrQuad(i,0);
        int nodeid2 = ElemArrQuad(i,1);
        int nodeid3 = ElemArrQuad(i,2);
        int nodeid4 = ElemArrQuad(i,3);
        int nodeid5 = ElemArrQuad(i,4);
        int nodeid6 = ElemArrQuad(i,5);
        int nodeid7 = ElemArrQuad(i,6);
        int nodeid8 = ElemArrQuad(i,7);
        Eigen::MatrixXd nodeIds(1, 8);
        nodeIds << nodeid1, nodeid2, nodeid3, nodeid4, nodeid5, nodeid6, nodeid7, nodeid8; //comma initialization
        //std::cout << "nodeIds" << nodeIds << std::endl;
       
        for (int m=0;m<NumNodesPerElem;m++)
        {   
            for (int n=0;n<NumNodesPerElem;n++)
            {
                int globalI1 = 3*ElemArrQuad(i,m);
                int globalI2 = 3*ElemArrQuad(i,m)+1;
                int globalI3 = 3*ElemArrQuad(i,m)+2;

                int globalJ1=3*ElemArrQuad(i,n);
                int globalJ2=3*ElemArrQuad(i,n)+1;
                int globalJ3=3*ElemArrQuad(i,n)+2;
               
                KGlobal(globalI1,globalJ1) = KGlobal(globalI1,globalJ1) + Klocal(3*(m),3*(n));
                KGlobal(globalI1,globalJ2) = KGlobal(globalI1,globalJ1) + Klocal(3*(m),3*(n)+1);
                KGlobal(globalI1,globalJ3) = KGlobal(globalI1,globalJ1) + Klocal(3*(m),3*(n)+2);
                KGlobal(globalI2,globalJ1) = KGlobal(globalI1,globalJ1) + Klocal(3*(m)+1,3*(n));
                KGlobal(globalI2,globalJ2) = KGlobal(globalI1,globalJ1) + Klocal(3*(m)+1,3*(n)+1);
                KGlobal(globalI2,globalJ3) = KGlobal(globalI1,globalJ1) + Klocal(3*(m)+1,3*(n)+2);
                KGlobal(globalI3,globalJ1) = KGlobal(globalI1,globalJ1) + Klocal(3*(m)+2,3*(n));
                KGlobal(globalI3,globalJ2) = KGlobal(globalI1,globalJ1) + Klocal(3*(m)+2,3*(n)+1);
                KGlobal(globalI3,globalJ3) = KGlobal(globalI1,globalJ1) + Klocal(3*(m)+2,3*(n)+2);


                //KGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) = KGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) + Klocal_arr[i](m,n);
                //MGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) = MGlobal(int(nodeIds(0,m)),int(nodeIds(0,n))) + Mlocal_arr[i](m,n);
            }
            FGlobal(int(nodeIds(0,m)),0) = FGlobal(int(nodeIds(0,m)),0) + Flocal(m,0);
        }
   
    }

    //Assembly
    for (int i=0;i<NumOfElems;i++)
    {
         }
   //std::cout << "KGlobal " << KGlobal << std::endl;    
    return 0;
}

Eigen::MatrixXf GetDGlobalArr(int NElemsX, int NElemsY, int NElemsZ, int order, int dim, int bondaryCondition, float F)
{
    
    int NumOfElems = NElemsX*NElemsY*NElemsZ;
    int NPerElems =8; //nodes per Elem

    float xMin = 0.0;
    float yMin = 0.0;
    float zMin = 0.0;
    float xMax = 0.1;
    float yMax = 0.1;
    float zMax = 1.0;
    double g2 = 0.05; //displacement at z=zMax
    int  NumNodesTotal  = (NElemsX+1)*(NElemsY+1)*(NElemsZ+1);

    if (order==1) 
    {
        NumNodesTotal = (NElemsX+1)*(NElemsY+1)*(NElemsZ+1);//linear
        NPerElems = 8;
    }
    
    //createmesh
    Eigen::MatrixXd  ElemArrQuad = GetElemArr(NElemsX, NElemsY, NElemsZ, NPerElems);
    Eigen::MatrixXf GNodeCordsArr = GetGNodeCordsArr(NElemsX, NElemsY,NElemsZ, NPerElems ,dim, xMin, yMin, zMin, xMax, yMax, zMax);
  
    Eigen::MatrixXf  KGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,dim*NumNodesTotal);
    Eigen::MatrixXf  MGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,dim*NumNodesTotal);
    Eigen::MatrixXf  FGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,1);
    Eigen::MatrixXf  DGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,1);

    std::cout << "GetKMandF3DFE\n";
    //Getting KGlobal, MGlobal, FGlobal
   GetKMandF3DFE(NElemsX, NElemsY, NElemsZ, GNodeCordsArr, ElemArrQuad, order, bondaryCondition, F, KGlobal, MGlobal, FGlobal);
    std::cout << "GetKMandF3DFEdone\n";
// Dirichlet BC
   // std::cout << "GNodeCordsArr," << GNodeCordsArr << std::endl;
    std::cout << "KGlobal," << KGlobal << std::endl;
    int nDirichletnodes = (NElemsX+1)*(NElemsY+1)*dim+(NElemsX+1)*(NElemsY+1);
    Eigen::MatrixXf  DirichBC = Eigen::MatrixXf::Zero(nDirichletnodes,2);
    //Eigen::MatrixXf  DirichBC;
    int j =0;
    for (int i=0;i<NumNodesTotal;i++)
    { 
        if(GNodeCordsArr(i,2)==0)   //at z=0
        {
            DirichBC(j,0)=3*i;
            DirichBC(j,1)=0;
            DirichBC(j+1,0)=3*i+1;
            DirichBC(j+1,1)=0;
            DirichBC(j+2,0)=3*i+2;
            DirichBC(j+2,1)=0;
            j=j+3;
        }
        if(GNodeCordsArr(i,2)==zMax)   //at z=zMax
        {
            DirichBC(j,0)=3*i+1;
            DirichBC(j,1)=g2;
            j=j+1;
        }
    }
    std::cout << "The vector DirichBC is of size " << DirichBC.size() << std::endl;
    //for 1 elem
    int nNotDirchiletnodes = dim*NumNodesTotal-nDirichletnodes;
    std::cout << "DirichBC " << DirichBC << std::endl;  
    //Applying Dirichilet Boundary conditions

//      //Apply Dirichlet's BCs
//    3D Hexahedral Mesh: ([0, 0.1]x[0, 0.1]x[0, 1.0], use a 4 x 4 x 40 element mesh.)
          //  Fully constrained on z = 0m face: g = 0 along z = 0m
        //    Y-displacement on z = 1m face: g2 = 0.05 along z = 1.0m
  
    //KDGlobal, MDGlobal and FDGlobal Dirichlet
    Eigen::MatrixXf  KDGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal-nDirichletnodes,dim*NumNodesTotal-nDirichletnodes);
   // Eigen::MatrixXf  MDGlobal = Eigen::MatrixXf::Zero(nNotDirchiletnodes,nNotDirchiletnodes);
    Eigen::MatrixXf  FDGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal-nDirichletnodes,1);
    Eigen::MatrixXf  NotDirchNodes = Eigen::MatrixXf::Zero(dim*NumNodesTotal-nDirichletnodes,1);
  std::cout << "DirichBC1 " << std::endl;  
  
    //Non Dirchlet Nodes matrix
    int count=0;
    int flag=0;
    for (int i=0;i<dim*NumNodesTotal;i++)
    {
        flag=0;
        for (int k=0;k<nDirichletnodes;k++)
        {
            if(i==DirichBC(k,0))
            {
                flag=1;
            }
        }

        if(flag==0)
        {
            NotDirchNodes(count,0)=i;
            count++;
        }
    }
    std::cout << "NotDirchNodes " << NotDirchNodes << std::endl;  
  
    //KDGlobal
    for (int k=0;k<nNotDirchiletnodes;k++)
    {
        for (int l=0;l<nNotDirchiletnodes;l++)
        {
            KDGlobal(k,l)=KGlobal(int(NotDirchNodes(k,0)),int(NotDirchNodes(l,0)));
        }
    }
    std::cout << "KDGlobal " << std::endl;  

    //fdash
    Eigen::MatrixXf  F_forSubstract = Eigen::MatrixXf::Zero(nNotDirchiletnodes,1);
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        for (int j=0;j<nDirichletnodes;j++)
        {
            F_forSubstract(i,0)=KGlobal(int(NotDirchNodes(i,0)),int(DirichBC(j,0)))*DirichBC(j,1) + F_forSubstract(i,0);
        }
    }
    std::cout << "F_forSubstract " << std::endl;  
 
    // FDGlobal
    for (int i=0;i<nNotDirchiletnodes;i++)
    {
        FDGlobal(i,0)=FGlobal(int(NotDirchNodes(i,0)),0) - F_forSubstract(i,0);
    }
   std::cout << "FDGlobal " << FDGlobal << std::endl;    
   std::cout << "DGlobaldone1 " << std::endl;    

    Eigen::MatrixXf  DGlobal_D = Eigen::MatrixXf::Zero(dim*NumNodesTotal-nDirichletnodes,1);
    DGlobal_D = KDGlobal.inverse() * FDGlobal;
   std::cout << "DGlobaldone2 " << std::endl;    
   std::cout << "DirichBC " << DirichBC << std::endl;    
   std::cout << "NotDirchNodes " << NotDirchNodes << std::endl;    
   std::cout << "DGlobal " << DGlobal << std::endl;    
   std::cout << "DGlobal " << DGlobal << std::endl;    

    for (int m=0;m<nDirichletnodes;m++)
    { 
        std::cout << "int(DirichBC(m,0) " << int(DirichBC(m,0)) << std::endl;    
        std::cout << "DGlobal(int(DirichBC(m,0)),0) " << DGlobal(int(DirichBC(m,0)),0) << std::endl;    
        std::cout << "DirichBC(m,1) " << DirichBC(m,1) << std::endl;    
        DGlobal(int(DirichBC(m,0)),0)=double(DirichBC(m,1));
        std::cout << "DGlobal(int(DirichBC(m,0)),0) " << DGlobal(int(DirichBC(m,0)),0) << std::endl;    
    }
   std::cout << "DGlobaldone3 " << std::endl;    

    for (int n=0;n<nNotDirchiletnodes;n++)
    {   
        DGlobal(int(NotDirchNodes(n,0)),0)=DGlobal_D(n,0);
    }

   std::cout << "DirichBC " << DirichBC << std::endl;    
   std::cout << "DirichBC " << DirichBC << std::endl;    
   std::cout << "DGlobal " << DGlobal << std::endl;    

    return DGlobal;
}


int main()
{
    //Inputs initialization
    double E = 1000 ;//Pa
    double Nu = 0.3 ;//Pa
    double g1 = 0;//, 
    double g2 = 0.05;// m, 

    int NElemsX =4;//horizontal
    int NElemsY =4;// vertical
    int NElemsZ =40;// Z direction
    NElemsX =1;//horizontal
    NElemsY =1;// vertical
    NElemsZ =1;// Z direction
    int NumOfElems = NElemsX*NElemsY*NElemsZ;
    int NumNodesTotal = (NElemsX+1)*(NElemsY+1)*(NElemsZ+1);
    int NPerElems =8;
    int dim = 3;
    float xMin = 0.0;
    float yMin = 0.0;
    float zMin = 0.0;
    float xMax = 0.1;
    float yMax = 0.1;
    float zMax = 1.0;


    float lambda = E/(2*(1+Nu));
    float Mu = E*Nu/((1+Nu)*(1-2*Nu));
    int order = 1;
    //Apply Dirichlet's BCs
    //     /* for the following sets of boundary conditions and forcing functions:
    //        1. (2D Quadrilateral Mesh): (x ∈ [0, 1], y ∈ [0, 1], use a 20 x 20 element mesh.)
    //        g = 300 K along x = 0 m (left edge) and g = 310 K along x = 1.0 m (right edge)
    //        Initial condition: u(x, y, 0) = 300K for x < 0.5m and u(x, y, 0) = 300 + 20(x − 0.5)K for x ≥ 0.5m

    //     where  f_tilde and f_bar are constant values.*/
    //     //κij = ¯κδij , where ¯κ = 385 watt.m−1K−1, ρ = 6000 kg.m−3 and ρc = 3.8151 x 106 N.m−2K−1

    std::cout << "write data to the file\n";
    std::ofstream out ("HW4ElemArr.txt");
    //createmesh for output
    Eigen::MatrixXd  ElemArrQuad = GetElemArr(NElemsX, NElemsY, NElemsZ, NPerElems);
    //out << "ElemArrQuad" "\n";
    std::cout << ElemArrQuad << std::endl;
    out << ElemArrQuad << "\n";
    out.close();

    std::ofstream out1 ("HW4NodeCoordsArr.txt");
    //create node coordinates array for output
    Eigen::MatrixXf GNodeCordsArr = GetGNodeCordsArr(NElemsX, NElemsY,NElemsZ, NPerElems ,dim, xMin, yMin, zMin, xMax, yMax, zMax);
    //out1 << "GNodeCordsArr" "\n";
    out1 << GNodeCordsArr << "\n";
    out1.close();


    Eigen::MatrixXf  KGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,dim*NumNodesTotal);
    Eigen::MatrixXf  MGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,dim*NumNodesTotal);
    Eigen::MatrixXf  FGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,1);
    Eigen::MatrixXf  DGlobal = Eigen::MatrixXf::Zero(dim*NumNodesTotal,1);

    float F=0.0;
    int bondaryCondition = 1;
    //Getting DGlobal
    std::cout << "DGlobal\n";
    
    std::ofstream out2 ("HW4DGlobalArr.txt");
    Eigen::MatrixXf DGlobalArr = GetDGlobalArr(NElemsX, NElemsY,NElemsZ, order, dim, bondaryCondition, F);
    //out2 << "DGlobalArr" "\n";
    out2 << DGlobalArr << "\n";
    out2.close();
    std::cout << "DGlobalDone\n";
    std::cout << ElemArrQuad << std::endl;
    return 0;
}
