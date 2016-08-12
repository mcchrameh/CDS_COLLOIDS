#include "CDS_BASE.h"
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

//#include "mpi.h"
#define MAX_Partices  100000
#define Random_min  -0.05
#define Random_max  0.05 
#define MAX_LINE_LENGTH 200
#define P2(i,j)    P2[i][j]
#define CP4(i,j)  CP4[i][j]
BASE::BASE()
     {
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);      
       double rationals[26];
       int integers[8];
       if(rank==0)
	  {
	    Read_input_parameters(integers, rationals);
	   }
       MPI_Bcast(integers,8,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(rationals,26,MPI_DOUBLE,0,MPI_COMM_WORLD);
       Nx   = integers[0];     Ny    = integers[1]; //Nz=integers[2];
       Procx=integers[2];      Procy = integers[3]; //Procz =integers[5];
       Nbig =integers[4];      Nsmall= integers[5];
       Nparticles=integers[6];
       M = rationals[0]; u=rationals[1]; B = rationals[2];
       D = rationals[3]; A=rationals[4]; v = rationals[5];
       tau= rationals[6]; F=rationals[7]; r=rationals[8];
       Max_time_iteration=rationals[9];
       U1=rationals[10]; CSI=rationals[11];
       ENNE=rationals[12];      ALPHA =rationals[13];
       beta =rationals[14];     GAMMs =rationals[15];
       GAMMb =rationals[16];    temp1 =rationals[17];
       R1small=rationals[18];   R1big =rationals[19];
       ALP1 =rationals[20];     dxx   =rationals[21];
       EMME =rationals[22];     sigma =rationals[23];
       P00s =rationals[24];     P00b  =rationals[25];

        delta_t=0.05; 
       //Nbig = integers[4];
       //Nparticles=integers[5];
       RR1small= R1small*pow((1.0 +(1.0/log(2.0))),(1.0/ALP1));
       RR1big=R1big*pow((1.0 +(1.0/log(2.0))),(1.0/ALP1));
       RCUT = RR1big*0.01; //<---This is arbitrary
       printf("RCUT=%lf\n",RCUT);

      int coords[2], dims[2], periods[2];

      MPI_Comm new_comm =CreatCartesianTopology();
      MPI_Comm_rank(new_comm,&my2drank);
      MPI_Cart_get(new_comm,2,dims,periods,coords);
      assert(Procx*Procy==size);
      printf("dimensions of topology=(%d,%d)\n",dims[0],dims[1]);
      
      if(Nx%Procx!=0)
	{
         if(coords[0]<(dims[0]-1))
           {
             nlocalx=(int)Nx/Procx;
	   }
	 else
          {
            nlocalx = (int)Nx/Procx + (int)Nx%Procx;

          }
	}
      else
	{
		
          nlocalx = (int)Nx/Procx;
        }

      if(Ny%Procy!=0)
	{
	  if(coords[1]<(dims[1]-1))
	    {
	      nlocaly = (int)Ny/Procy;
	    }
	 else
	    {
	      nlocaly = (int)Ny/Procy + (int)Ny%Procy;
	    }
	 }
      else
	 {
             nlocaly = (int)Ny/Procy;

          }
   



/*--Compute Number of linked cells in x and y direction -------*/
      LinkedCells.resize(2);
      CellSize.resize(2);
      for (int direction=0; direction<2; direction++) 
          {
	    LinkedCells[direction] = (int) (al[direction]/RCUT) ; 
	    CellSize[direction]    = (double)(al[direction]/LinkedCells[direction]);
	    printf("Cellsize[%d]=%lf\n",direction, LinkedCells[direction]);
          }
      /*---Add ghost cells to number of linked cells---*/
      for (int direction=0; direction<2; direction++)
	   {
              LinkedCells[direction]=LinkedCells[direction] + 2;
	   }

      head.resize(LinkedCells[0]*LinkedCells[1]);
      LinkedCell_list.resize(Nparticles);
    //-----------------------------------------------------------------
    //Establish the local number of particles on each MPI process
    //---------------------------------------------------------------
           /*
	   nbig_local= (int)(Nbig/size);
	   nbig_local=(my2drank<Nbig%size)?nbig_local +1 : nbig_local ; //<-----check this
	   nsmall_local= (int)(Nsmall/size);
	   nsmall_local=(my2drank<Nsmall%size)?nsmall_local +1 : nsmall_local ;

           nlocal_particles =Nparticles/size;
           N_start = (int)my2drank*nlocal_particles;
	   if(Nparticles!=(Nbig + Nsmall))
	   {
            printf("particles don't match; Nbig=%d, Nsmall=%d, Nparticles=%d\n", Nbig, Nsmall, Nparticles);
	    printf("Condition Nparticles = Nbig + Nsmall, is not satisfied");
	    assert(Nparticles==(Nbig + Nsmall));
	   }
	   */
           nlocal_particles =Nparticles/size;
           N_start = (int)my2drank*nlocal_particles;

           

	   if(Nparticles%size)
	    {
             nlocal_particles=(my2drank<Nparticles%size)?nlocal_particles +1 : nlocal_particles ;
	     N_start =N_start + ((my2drank < ((Nbig + Nsmall)%size) ) ? my2drank : ((Nbig + Nsmall)%size));
	    }
         	    
            


  recvcounts=new int[size];
  displacement=new int[size]; 
  for (int i = 0; i < size; i++ )
     { 
        recvcounts[i] = (( i < (Nparticles%size) )? nlocal_particles + 1: nlocal_particles)*sizeof(BODY)/sizeof(double);
        displacement[i] = ((int)(i*nlocal_particles) + ((i < Nparticles%size ) ? i : (Nparticles%size)))*sizeof(BODY)/sizeof(double);
         //printf("%d, %d, %d, %d\n", rank, i, recvcounts[i], displacement[i]);
     }
 
	   
//------------------------------------------------------------------------------------

 

       PHI        = new double*[nlocalx+2];
       PHI_old    = new double*[nlocalx+2];
       gamma      = new double*[nlocalx+2];
       Laplacian2 = new double*[nlocalx+2];
       PP3        = new double*[nlocalx+2];
       P2         = new double*[nlocalx+2];
       CP4        = new double*[nlocalx+2];
       PPP        = new double*[nlocalx+2];
       P3         = new double*[nlocalx+2];
       PHI_local_result = new double[nlocalx*nlocaly];
       matrix_lower = new double[nlocaly];
       matrix_upper = new double[nlocaly];
       matrix_left  = new double[nlocalx+2];
       matrix_right = new double[nlocalx+2];
      // std::shared_ptr<Node> nd( new Node);//malloc(sizeof(Node));

       for(int i=0;i<nlocalx+2; i++)
	  {
	   PHI[i]=new double[nlocalx+2];
	   PHI_old[i]=new double[nlocaly+2]; 
	   gamma[i]=new double[nlocaly+2];
	   Laplacian2[i]=new double [nlocaly+2];
	   PP3[i]  =new double [nlocaly+2];
	   P2[i]   = new double [nlocaly+2];
	   CP4[i]  = new double [nlocaly+2];
	   PPP[i]  = new double [nlocaly+2];
	   P3[i]   = new double [nlocaly+2];


          }
      

       for(int i=0;i<nlocalx +2; i++)
	   {
            for(int j=0;j<nlocaly+2;j++)
	    {
               PHI[i][j]=0.0;
	       PHI_old[i][j]=0.0;
	       CP4[i][j]=0.0;
	       PP3[i][j]=0.0;
	       P2[i][j]=0.0;
	       gamma[i][j]=0.0;
	       Laplacian2[i][j]=0.0;
	       PPP[i][j]=0.0;
	       P3[i][j]=0.0;

	    }

	   }


//   XI.resize(nlocal_particles);
   bd.resize(Nparticles);   //= new BODY[Nparticles];
   //Recv_bodies= new BODY[nlocal_particles];
   //Send_bodies= new BODY[nlocal_particles];
  // bd.resize(Nparticles);i
   send_bodies_left.resize(nlocal_particles);
   send_bodies_right.resize(nlocal_particles);
   send_bodies_upper.resize(nlocal_particles);
   send_bodies_lower.resize(nlocal_particles);


   Recv_bodies.resize(nlocal_particles);
   send_bodies.resize(nlocal_particles);
   R12.resize(Nparticles);
   Dx.resize(Nparticles);
   Dy.resize(Nparticles);
   RI.resize(Nparticles);
   Dx1.resize(Nx);
   Dy1.resize(Nx);
   //PPP.resize(nlocalx);
   POS.resize(Nx);//<--- rethink about the size of this, it is bether to use a linked list
   VerletList_local.resize(Nparticles);
   VerletTable_local.resize(Nparticles);
   TOTX.resize(Nparticles);
   TOTY.resize(Nparticles);
   for(int i=0;i<Nx;i++)
      {
	POS[i].resize(Ny);
	Dx1[i].resize(Ny);
	Dy1[i].resize(Ny);
      }

   for(int i=0; i<Nparticles;i++)
      {
        Dx[i].resize(Nparticles);
        Dy[i].resize(Nparticles);
        R12[i].resize(Nparticles);
        RI[i].resize(Nparticles);
        VerletTable_local[i].resize(Nparticles);

      }
   
    blocklengths[0] = 1;//{1,1,1,1,1,1};
    blocklengths[1] = 1;
    blocklengths[2] = 1;
    blocklengths[3] = 1;
    blocklengths[4] = 1;
    blocklengths[5] = 1;
    types[0] = MPI_INT;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_DOUBLE;
    types[4] = MPI_DOUBLE;
    types[5] = MPI_DOUBLE;
    //{MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    offsets[0] = offsetof(BODY, index);
    offsets[1] = offsetof(BODY, r[0]);
    offsets[2] = offsetof(BODY, r[1]);
    offsets[3] = offsetof(BODY, v[0]);
    offsets[4] = offsetof(BODY, v[1]);
    offsets[5] = offsetof(BODY, mass);// {offsetof(BODY, index),offsetof(BODY, x),offsetof(BODY, y),offsetof(BODY, vx),offsetof(BODY, vy),                  offsetof(BODY, mass)};
        //MPI_DOUBLE   
      //MPI_Type_create_struct(6,blocklengths,offsets,types,&mpi_body_type);
    MPI_Type_create_struct(6,blocklengths,offsets,types,&mpi_body_type);
    MPI_Type_commit(&mpi_body_type);
 
  printf("Base Construction made. \n");
  FreeCommunicator(new_comm); 

 }

BASE::~BASE()
     {
        for(int i=0;i<nlocalx+2;i++)
           {
             delete [] PHI[i];
	     delete [] PHI_old[i];
	     delete [] gamma[i];   
	     delete [] Laplacian2[i];
	     delete [] PP3[i];
	     delete [] P2[i];
	     delete [] CP4[i];
	     delete [] PPP[i];
	     delete [] P3[i];

	   }
	//for(int i=0;i<Nx;i++)
	  // {
           //  delete [] PP3[i];

	  // }
	  delete [] PHI;
	  delete [] PHI_old;     
	  delete [] Laplacian2;
	  delete [] gamma;	    
	  delete [] P3;
	 // delete [] bd;
	 // delete [] Recv_bodies;
	 // delete [] Send_bodies;
	  delete [] recvcounts;
	  delete [] displacement;
	  delete [] PP3;
	  delete [] P2;
	  delete [] CP4;
	  delete [] PPP;
	  delete [] PHI_local_result;

	  delete [] matrix_lower;
	  delete [] matrix_upper;
	  delete [] matrix_left;
	  delete [] matrix_right;
	 // delete  nd;
	  std::cout<<"Destructor in base class called"<<std::endl;
	//  std::cout<<"Destructor in base class called: colloids destroyed"<<std::endl; 



     }
void BASE::InitializePolymer(MPI_Comm new_comm)
    {
      MPI_Comm_rank(new_comm, &my2drank);
      MPI_Comm_size(new_comm, &size);
	        
       srand(time(NULL) + my2drank);

        for(int i=0; i<nlocalx;i++)
           {
	    for(int j=0;j<nlocaly;j++)
	       {
	         double range =Random_max-Random_min;
	         double div =RAND_MAX / range;
	         PHI_old[i][j]=Random_min + (rand()/div);
	        }
	    }

     }





void BASE::Read_input_parameters(int *integers, double *rationals)
     {
       //int value(0);
       FILE* file;
       char Data[MAX_LINE_LENGTH], *string;
	    if((file=fopen("ParameterFile.dat","r"))==NULL)
              {
	         printf("Error opening ParameterFile.dat\n");
	         return;
	      }
        
	string=fgets(Data,MAX_LINE_LENGTH,file);
          		
	int value=fscanf(file,"%d\n",&integers[0]);

	 string=fgets(Data,MAX_LINE_LENGTH,file);
         
	 value=fscanf(file,"%d\n",&integers[1]);
	 
         string=fgets(Data,MAX_LINE_LENGTH,file);
	 
	 value=fscanf(file,"%d\n",&integers[2]);
	 
	 string=fgets(Data,MAX_LINE_LENGTH,file);
         
	 value=fscanf(file,"%d\n",&integers[3]);
	 
	 string=fgets(Data,MAX_LINE_LENGTH,file);
	 
	 value=fscanf(file,"%d\n",&integers[4]);
	 
	string=fgets(Data,MAX_LINE_LENGTH,file);
         
	value=fscanf(file,"%d\n",&integers[5]);
	 
	string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%d\n",&integers[6]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%d\n",&integers[7]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
 	value=fscanf(file,"%lf\n",&rationals[0]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[1]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[2]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[3]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[4]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[5]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[6]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[7]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[8]);
        string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%lf\n",&rationals[9]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[10]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
	value=fscanf(file,"%lf\n",&rationals[11]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[12]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[13]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[14]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[15]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[16]);

	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[17]);

	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[18]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[19]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[20]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[21]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[22]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[23]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[24]);
	string=fgets(Data,MAX_LINE_LENGTH,file);
        value=fscanf(file,"%lf\n",&rationals[25]);
//	string=fgets(Data,MAX_LINE_LENGTH,file);
//        value=fscanf(file,"%lf\n",&rationals[26]);
	
	
       // fgets(Data,MAX_LINE_LENGTH,file);
	//printf("value=%d,character=%s\n",value,string[0]);
	if((value!=0)||(string!=NULL))
	  {
           printf("failed to read in\n");
	  }
        fclose(file);
 
}



MPI_Comm BASE::CreatCartesianTopology()
      {
         MPI_Comm   new_comm;
         if(Procx*Procy!=size)
           {
             printf("Number of processors is not factorizable. That is Px*Py=%d, is not equal to the number of processors=%d\n ",Procx*Procy,size);
             printf("Hint: Check Parameter.dat file to change the number of cores\n");

		
	    exit(0);
	  }
        const int ROWS=0, COL=1;	     
        int periods[2];
        periods[0]=1; periods[1]=1;
        int coords[2];
        int dims[2];
        dims[ROWS]=Procx;
        dims[COL] =Procy;
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Cart_create( MPI_COMM_WORLD, 2,dims, periods, 0, &new_comm );
        MPI_Comm_rank(new_comm,&my2drank);
        MPI_Cart_coords(new_comm,my2drank,2,coords);
        MPI_Comm_size( new_comm, &size );

         return new_comm;

 }

void BASE::GenerateColloids(MPI_Comm new_comm)
{
   

   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   double c[2],gap[2];
  // double origAtom[3][2] = {{0.0,  0.0},
  //                       { 0.0, 0.5}, { 0.5, 0.0}};//the number of FCC atoms in a cell in 2D is 3 ie 4*0.5 + 4*.25. 4 atoms at the corner, each shared by 4 neighbors and 4 atoms on the edge each shared by two neigbors
   srand(time(NULL) + my2drank);
   
    for(int a=0; a<2; a++)
       {
	 gap[a] = al[a]/(double)InitUcell[a];
       }
     int i=0;

            for( i=0; i<nlocal_particles; i++) 
	       {
		//for(a=0; a<2; a++)
	          double range =1.0/4.0;
                    double div =RAND_MAX / range;
                 bd[i].r[0]=  ((double)rand()/(div))+0.01;//gap[0]*origAtom[j][0]
                 bd[i].r[1]=  ((double)rand()/(div))+0.01;//gap[1]*origAtom[j][1];
                 bd[i].v[0]=0.0;
                 bd[i].v[1]=0.0;
                 bd[i].index=1;
                 printf("Colloid (rank=%d, dx=%lf, dy=%lf) \n",my2drank, bd[i].r[0], bd[i].r[1]);
                
               }
	 
      
 /* Total # of atoms summed over processors */
 int nglob=0;
 MPI_Allreduce(&i,&nglob,1,MPI_INT,MPI_SUM,new_comm);
 if (my2drank == 0) printf("nglob = %d\n",nglob);
}


int BASE::OnBoundary( int ku,int i)
{
   int kd,kdd;
   kd = ku/2; /* x(0)|y(1)|z(2) direction */ 
   kdd = ku%2; /* Lower(0)|higher(1) direction */
   if (kdd == 0)
      return bd[i].r[kd]< RCUT;
   else if (kdd==1)
      return al[kd]-RCUT < bd[i].r[kd];
   else
      return 0;
	  
}




void BASE::SendBodies( MPI_Comm new_comm)
{

   MPI_Status status;
   MPI_Request request;	   
   /* ku represents the four corners/neighbors of a process */
   int ku,ku_0,ku_1,nsd;
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   int N_recv=0, tag=1;
   for(int kd=0;kd<2;kd++) //loop over x and y direction
      {
       for (int bc=0;bc<2;bc++)  //bc=0 for lower(left) boundary, bc=1 for upper(right) boundary 
	    {
              	    
	      lsb[2*kd+bc][0]=0; //reset the number of to be copied atoms on all 4 boundaries
            }
       /* scan all atoms both resident and copies to identify boundary atoms */
       for(int i=0;i<nlocal_particles+N_recv;i++)
	   {
             //MPI_Cart_shift(new_comm, kd, 1,ku_0, ku_1);
             for(int bc=0;bc<2;bc++)
		{
	          ku=2*kd + bc;
	         if (OnBoundary(ku,i))
		    {
                      lsb[ku][++(lsb[ku][0])]=i;
		     // sendBodies.push_back(bd[i]);
		    }
	          }
          
             }
      }
         printf("Number of Boundary atoms,left=%d\n",lsb[0][0] );
        //left boundary buffer
         
	 send_bodies_left.resize(lsb[0][0]);
	 for(int i=0;i<lsb[0][0];i++)
            {
	     send_bodies_left.push_back(bd[i]);
	    }
	 //right boundary buffer
	 send_bodies_right.resize(lsb[2][0]);
	 for(int i=0;i<lsb[2][0];i++)
	    {       
	     send_bodies_right.push_back(bd[i]);
            }
	 //upper boundary buffer
	 send_bodies_upper.resize(lsb[1][0]);
	 for(int i=0;i<lsb[1][0];i++)
	    {       
               send_bodies_upper.push_back(bd[i]);
            }
	 //lower boundary buffer
	 send_bodies_lower.resize(lsb[3][0]);
	 for(int i=0;i<lsb[3][0];i++)
	    {       
	      send_bodies_lower.push_back(bd[i]);
            }

           //get neighbors in x and y direction to send and recieve
	    MPI_Comm_rank(new_comm,&my2drank);
	    MPI_Cart_shift(new_comm,0, 1,&ku_0, &ku_1);
	   // MPI_Sendrecv_replace(&send_bodies_upper, lsb[1][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_1, tag,
	   //                       ku_0,tag, new_comm,&status);
	    MPI_Isend(&send_bodies_upper,lsb[1][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
	    MPI_Irecv(&send_bodies_upper,lsb[1][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);
            int new_size=bd.size() + lsb[1][0];
	    bd.resize(new_size);
	    bd.insert(std::end(bd),std::begin(send_bodies_upper),std::end(send_bodies_upper) );

	    //MPI_Sendrecv_replace(&send_bodies_lower, lsb[3][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_0, tag,
	   //	                    ku_1,tag, new_comm,&status);

	    MPI_Isend(&send_bodies_lower,lsb[3][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);
            MPI_Irecv(&send_bodies_lower,lsb[3][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
            int new_size1=bd.size() + lsb[3][0];
	    
	    bd.resize(new_size1);
	    bd.insert(std::end(bd),std::begin(send_bodies_lower),std::end(send_bodies_lower) );

	    MPI_Comm_rank(new_comm,&my2drank);
	    MPI_Cart_shift(new_comm,1, 1,&ku_0, &ku_1);
	    //MPI_Sendrecv_replace(&send_bodies_left, lsb[0][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_0, tag,
	    //                         ku_1,tag, new_comm,&status);
	    MPI_Isend(&send_bodies_left,lsb[0][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);
	    MPI_Irecv(&send_bodies_left,lsb[0][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
            int new_size2=bd.size() + lsb[0][0];   
            bd.resize(new_size2);
	    bd.insert(std::end(bd),std::begin(send_bodies_left),std::end(send_bodies_left) );
	                
	   // MPI_Sendrecv_replace(&send_bodies_right, lsb[2][0]*sizeof(BODY)/sizeof(double),  MPI_DOUBLE, ku_1, tag,
	    //                        ku_0,tag, new_comm,&status);
	    MPI_Isend(&send_bodies_right,lsb[2][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_1,tag,new_comm,&request);
	    MPI_Irecv(&send_bodies_right,lsb[2][0]*sizeof(BODY)/sizeof(double),MPI_DOUBLE,ku_0,tag,new_comm,&request);

	    int new_size3=bd.size() + lsb[2][0];   
            bd.resize(new_size3);
	    bd.insert(std::end(bd),std::begin(send_bodies_right),std::end(send_bodies_right) );

      
}


void BASE::SetEffectiveRadius(MPI_Comm new_comm)
{
  // printf("IN EFFECTIVE RAD FUNCTION");
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
  
  for(int i=N_start;i<N_start + nlocal_particles-1;i++)
 //  for(int i=0;i<Nparticles;i++)
      {
	for(int j=i+1; j<N_start + nlocal_particles;j++) //Nparticles becos of verlet list
	   {
	    if((bd[i].index==0)&&(bd[j].index==0))
	      {
	       R12[i][j]=RR1small*2.0;

	      }
	    if((bd[i].index==1)&&(bd[j].index==1))
	      {
	       R12[i][j]=RR1big*2.0;
	      }
	    if((bd[i].index==1)&&(bd[j].index==0))
	    {
	      R12[i][j]=RR1big + RR1small;
	    }
	    if((bd[i].index==0)&&(bd[j].index==1))
	    {
	      R12[i][j]=RR1big + RR1small;

	    }
	 }
     }

}
void BASE::DistributeColloids(MPI_Comm new_comm)
{
   //MPI_Status status;
  
  MPI_Comm_rank(new_comm, &my2drank);
  MPI_Comm_size(new_comm, &size);
 //use cyclic distribution
//MPI_Allgatherv(bd + N_start, nlocal_particles*sizeof(BODY)/sizeof(double),MPI_DOUBLE, bd, recvcounts, displacement,MPI_DOUBLE, new_comm );
}

void BASE::GenerateVerlet(MPI_Comm new_comm)
{
    //check this implementation. Check if the total number of particles per process is Nparticles. 
//Nparticles becos, after the call of MPI_AllgatherV, all processes should have Nparticles and generate only a local verlet list.
printf("IN VERLET FUNCTION");
 MPI_Comm_rank(new_comm, &my2drank);

 MPI_Comm_size(new_comm, &size);
    for(int i=N_start; i<N_start+nlocal_particles;i++)
       {
         VerletList_local[i]=0;
         for(int j=i+1;j<Nparticles;j++)
            {
             double dx=bd[i].r[0]-bd[j].r[0];
             double dy=bd[i].r[1]-bd[j].r[1];
             dx=dx-nlocalx*(int((dx/(nlocalx/2)+3)/2)-1);
             dy=dy-nlocaly*(int((dy/(nlocaly/2)+3)/2)-1);
             double R2=sqrt(pow(dx,2)+pow(dy,2));
             if(R2<R12[i][j])
               {
                VerletList_local[i]=VerletList_local[i] +1;
                VerletTable_local[i][VerletList_local[i]]=j;
	        Dx[i][j]=dx; //distance between body i and body j
	 //       Dx[j][i]=-dx;
	        Dy[i][j]=dy;
	   //     Dy[j][i]=-dy;
	        RI[i][j]=R2;
	     //   RI[j][i]=R2;
	      }
	    }
      }
	

}
void BASE:: FreeCommunicator(MPI_Comm new_comm)
{
         MPI_Comm_free(&new_comm);
		    

}
void BASE::Coupling(MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   int coords[2];
   MPI_Cart_coords(new_comm,my2drank,2,coords);
  // std::vector<auto> integration_index; //<--stores indices for integration


   int N1(0), GT1(0), GT2(0), X, Y;
   double RX1(0.0),RY1(0.0), R2(0.0), B1(0.0),B2(0.0), D1(0.0), D2(0.0);
  for(int i=N_start;i<N_start+nlocal_particles;i++)
   //for(int i=0; i<nlocal_particles; i++)
      {
       if(bd[i].index==0)//for small colloids
	 {
          P00=P00s;
          N1=(int)RR1small+1;
	  TOTX[i]=0;TOTY[i]=0;
	  for( int j=-N1;j<N1;j++)
	     {
               for(int k=-N1;k<N1;k++)
		  {
	           GT1=(int)bd[i].r[0] + j;
		   GT2=(int)bd[i].r[1] + k;
		   if(GT1>nlocalx) // bc condition on a process in x direction
		     {
                       GT1=GT1-nlocalx*((int)(GT1/nlocalx+1)-1);
		     }
		   if(GT1<0)
		     {
                      GT1=GT1 + nlocalx;
		     }
		   X=std::abs((int)GT1);
                   X = (X-my2drank)/size;//,<--global to local index

		   if(GT2>nlocaly)
	             {
                      GT2=GT2-nlocaly*((int)(GT2/nlocaly + 1)-1);
		     }
		   if(GT2<0)
		     {
		      GT2=GT2+nlocaly;
		     }
		   Y=std::abs((int)GT2);
                   Y = (Y-my2drank)/size;                   
		   RX1=bd[i].r[0]-(double)GT1;
		   RX1=RX1-nlocalx*((int)((RX1/(nlocalx/2) + 3)/2)-1);//<----check this
		   RY1=bd[i].r[1]-(double)GT2;
		   RY1=RY1-nlocaly*((int)((RY1/(nlocaly/2) + 3)/2)-1);
		   R2=sqrt(pow(RX1,2)+pow(RY1,2));
		   if(R2<RR1small)
		   {
                    Dx1[X][Y]=RX1;
		    Dy1[X][Y]=RY1;
		    B1 = 1.0-pow(R2/RR1small, ALP1);
		    B2 =exp(1.0-(1.0/B1));
		    D1 =ALP1/pow(RR1small,ALP1);
		    D2 =-(1.0)*(D1*pow(R2, ALP1-2.0)/(pow(B1,2)))*B2;//<--check this
		    POS[X][Y]=i;
                    PP3[X][Y]=PP3[X][Y] + B2*(PHI_old[X][Y]-P00);
		    if(R2<R1small)
		      {
                        PPP[X][Y]=1.0;
		      }
		    TOTX[i]+= -sigma*RX1*D2*pow((PHI_old[X][Y]-P00),2);
		    TOTY[i]+= -sigma*RY1*D2*pow((PHI_old[X][Y]-P00),2);
                   }
        	  }       
               }
   
	    }
      if (bd[i].index==1) //big particles
       {
	 /*here, we should consider only particles that are in the current process. We start by checking
	  weather the particles are within the MPI process.  */
	  N1=0;
          P00=P00b;
          N1=(int)RR1big+1;
	  TOTX[i]=0.0;
	  TOTY[i]=0.0;
           for(auto j=-N1;j<N1;j++)
	      {
               for(auto k=-N1;k<N1;k++)
  		  {
                    GT1=(int)bd[i].r[0] + j;
		    GT2=(int)bd[i].r[1] + k;
		    if(GT1>(nlocalx)) // bc condition on a process in x direction
                     {
  		       break;
     		       //GT1=GT1-nlocalx*((int)(GT1/nlocalx+1)-1);
		      }
		    if(GT1<0)
		      {
		       break;	      
		       //GT1=GT1 + nlocalx;
		      }
		     X=std::abs((int)GT1);
		      X = (X-my2drank)/size;//,<--global to local index
     
		    if(GT2>(nlocaly))
		      {
		       break;	      
		       //GT2=GT2-nlocaly*((int)(GT2/nlocaly + 1)-1);
		      }
		    if(GT2<0)
		      {
		       break;
		       //GT2=GT2+nlocaly;
		      }
		    Y=std::abs((int)GT2);
                    Y = (Y-my2drank)/size;

		    RX1=bd[i].r[0]-(double)GT1;
		    if((RX1>(double)nlocalx)||(RX1<0.0) )
		      {
		       break;	      
		      // RX1=RX1-nlocalx*((int)((RX1/(nlocalx/2) + 3)/2)-1);//<----check this
		      }

		    RY1=bd[i].r[1]-(double)GT2;
		    if((RY1>(double)nlocaly)||RY1<0.0)
	            {
		     break;
	   	    // RY1=RY1-nlocaly*((int)((RY1/(nlocaly/2) + 3)/2)-1);
		    }
		    R2=sqrt(pow(RX1,2)+pow(RY1,2));
		    if(R2<RR1big)
	              {
	               Dx1[X][Y]=RX1;
	               Dy1[X][Y]=RY1;
	               B1 = 1.0-pow(R2/RR1big, ALP1);
	               B2 =exp(1.0-(1.0/B1));
	               D1 =ALP1/pow(RR1big,ALP1);
	               D2 =-(1.0)*(D1*pow(R2,ALP1-2.0)/(pow(B1,2)))*B2;//<--check this
		       //printf("X=%d, Y=%d\n", X,Y);
		       POS[X][Y]=i;
		      // printf("PP3[%d][%d]=%lf, P00=%lf\n",X,Y, PP3[X][Y],P00);
	               PP3[X][Y]=PP3[X][Y] + B2*(P3[X][Y]-P00);
		   //    printf("PP3[%d][%d]=%lf, P00=%lf\n",X,Y, PP3[X][Y],P00);
		       if(R2<R1big)
			 {
                           PPP[X][Y]=1.0;
			 }
	               TOTX[i]+= -sigma*RX1*D2*pow((P3[X][Y]-P00),2);
	               TOTY[i]+= -sigma*RY1*D2*pow((P3[X][Y]-P00),2);
		      // printf("Colloid (X =%d, Y=%d,P3[%d][%d]=%lf) \n",X, Y,X,Y, P3[X][Y]);
		    //   printf("Colloid (X =%d, Y=%d,TOTX=%lf,TOTY=%lf) \n",X, Y, TOTX[i],TOTY[i]);

	              }
       		   }
                }
              }
           } //end of particle loop 
   //computation equation 29

     
   for(int i=0;i<nlocalx;i++)
      {
       for(int j=0;j<nlocaly;j++)
	  {
           if(POS[i][j]!=0)
            {
             int z= POS[i][j];
	      if(bd[z].index==0)
		 {
                  P00=P00s;
		 }
	      else if(bd[z].index==1)
		 {
                   P00=P00b;
		 }
	      else
		 {
	          continue;
		 }
	    }
           P2[i][j]=sigma*2.0*PP3[i][j];
	  // printf("P2[%d][%d]=%lf\n", i,j, P2[i][j]);
	  }
      }
   //calculate laplacian for use with polymer
   /*
   double AP3(0.0), BP3(0.0);

   for(int i=1;i<=nlocalx;i++)
      {
       for(int j=1;j<=nlocaly;j++)
	  {
            AP3=(1.0/6.0)*(P2(i+1,j)+P2(i-1,j) + P2(i,j+1) + P2(i,j-1));
	    BP3=(1.0/12.0)*(P2(i-1,j-1) + P2(i-1,j+1) + P2(i+1,j-1) + P2(i+1,j+1));
	    CP4(i,j)=AP3 + BP3;
	  }

      }
      */

}
