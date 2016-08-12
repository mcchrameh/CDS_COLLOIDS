#ifndef _BASECLASS_
#define _BASECLASS_

//#include "CDS_Polymer.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <assert.h>
#include "mpi.h"
#include <memory>
#include <cstddef>
#define EMPTY -1
#define NBMAX 100000  //Maximum # of copied boundary atoms per neighbor.

typedef struct{                  
 int index; //determines if body is small or larg ie index=1 corresponds to large colloid and index=0 otherwise 
 double r[2];    
 //double y;     
 double v[2];
 //double vy; 
 double mass;
 }BODY;

typedef struct Node{
   double cx;
   double cy;
   double cmx;
   double cmy;
   double nmass;
   double xmin, xmax;
   double ymin, ymax;
   struct Node* nw;
   struct Node* ne;
   struct Node* sw;
   struct Node* se;
   BODY* bd;
   int has_childs;
 } Node;


//U1 energy scale, kbT
class BASE
{

  protected:
	    int Nx, Ny, Procx, Procy,nlocalx,nlocaly,my2drank,rank,size, Nbig,Nsmall, Nparticles, nlocal_particles, N_start;
	    int nbig_local,nsmall_local,  *recvcounts, *displacement;
	    double **PHI, **PHI_old, **Laplacian2, **gamma,**PP3,**P2,**CP4,**P3,  *PHI_local_result, **PPP;
            int InitUcell[2]={16,16}; // 
	    double al[2]={1.0,1.0};  //<---box length in x and y direction
	    double  RCUT;

	    int lsb[6][NBMAX];
	    std::vector<int> LinkedCells;  //<--number of linked cells in x and y direction
	    std::vector<double> CellSize;
	    std::vector<int> head;
	    std::vector<int> LinkedCell_list;
	    //BODY *bd;
	    std::vector<BODY> bd;
	    std::vector<BODY> Recv_bodies;
	    std::vector<BODY> send_bodies;
	    std::vector<BODY> send_bodies_left;
	    std::vector<BODY> send_bodies_right;
	    std::vector<BODY> send_bodies_upper;
	    std::vector<BODY> send_bodies_lower;



	    Node *nd;
	   // std::vector<BODY> bd;
	    double ALP1, RR1small, R1small, R1big, RR1big, P00, P00s,sigma, P00b,ENNE, ALPHA, GAMMs, GAMMb, U1, beta, delta_t, temp1, dxx, EMME, CSI ;
            int MPI_return;	    
	    double D, M, u, A , v, r, F, tau, Max_time_iteration, B;
	    std::vector<std::vector<double> > R12; //effective radius 
	    std::vector<int >VerletList_local;
	    std::vector<std::vector<int> >VerletTable_local;
	    std::vector<std::vector<double> > Dx;
	    std::vector<std::vector<double> > Dx1;
	    std::vector<std::vector<double> > Dy;
            std::vector<std::vector<double> > Dy1;
	    std::vector<std::vector<double> > RI;
	    //std::vector<std::vector<double> > r;
	  //  std::vector<std::vector<double> > PPP;
	    std::vector<std::vector<int>  > POS;
	    std::vector<double> TOTX;
	    std::vector<double> TOTY;
	    double *matrix_lower;
	    double *matrix_upper;
	    double *matrix_left;
	    double *matrix_right;
	    
	    const int nitems=5;
	    int  blocklengths[6] ;//= {1,1,1,1,1,1};
	     MPI_Datatype types[6] ;//= {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
             MPI_Datatype mpi_body_type;
             MPI_Aint     offsets[6];//={offsetof(BODY, index),offsetof(BODY, x),offsetof(BODY, y),offsetof(BODY, vx),offsetof(BODY, vy),offsetof(BODY, mass)};
	     //MPI_DOUBLE   
	    
	     
             //MPI_Type_create_struct(6,blocklengths,offsets,types,&mpi_body_type);
	    // MPI_Type_create_stuct(6,blocklengths,offsets,types,&mpi_body_type);
	    // MPI_Type_commit(&mpi_body_type);
            

  public:
	 BASE();
virtual ~BASE();
	void InitializePolymer(MPI_Comm new_comm);
	void Read_input_parameters(int *integers, double *rationals);
	MPI_Comm CreatCartesianTopology();
	void GenerateColloids(MPI_Comm new_comm);
	void DistributeColloids(MPI_Comm new_comm);
	void GenerateVerlet(MPI_Comm new_comm);
	void SendBodies( MPI_Comm new_comm);
        void SetEffectiveRadius(MPI_Comm new_comm);
	void  FreeCommunicator(MPI_Comm new_comm);
        void Coupling(MPI_Comm new_comm);
	int  OnBoundary( int ku,int i);
       // Node *CreateNode(double xim, double xmax,double ymin, double ymax);
       //void AddBodyToTree(Node* root, BODY* bd);
//       friend void Polymer2D::ExchangeData(MPI_Comm new_comm, double **array);		 






 };
#endif
