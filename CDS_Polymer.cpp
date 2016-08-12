#include "CDS_Polymer.h"
#include <cmath>



#define gamma(i,j)     gamma[i][j]
#define PHI_old(i,j)   PHI_old[i][j]
#define PHI(i,j)       PHI[i][j]
#define Laplacian2(i,j) Laplacian2[i][j]
#define CP4(i,j)  CP4[i][j]
#define P2(i,j)    P2[i][j]

Polymer2D::Polymer2D():BASE()
  { 
   


  }
Polymer2D::~Polymer2D()
 {
   std::cout<<"Polymer destructor called"<<std::endl;
 }
void Polymer2D::setLaplacianBase(MPI_Comm new_comm)
 {
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);    
   //cout<<"gamma(0,0)="<<gamma(0,0)<<"\t";
   double AP=0.0, BP=0.0,  ATP=0.0;//gtemp=0.0;
   for(int i=1;i<=nlocalx;i++)
      {
       for(int j=1;j<=nlocaly;j++)
          {
           AP=(1.0/6.0)*(PHI_old(i+1,j) + PHI_old(i-1,j)
                + PHI_old(i,j+1) + PHI_old(i,j-1));
           BP=(1.0/12.0)*(PHI_old(i-1,j+1) + PHI_old(i-1,j-1)
                + PHI_old(i+1,j+1) + PHI_old(i+1,j-1));
           ATP = AP + BP;
           gamma(i,j)=g(PHI_old(i,j))+ D*(ATP-PHI_old(i,j))-PHI_old(i,j);
           }
     }

    //for polymer colloid coupling
    double AP3(0.0), BP3(0.0);
    for(int i=1;i<=nlocalx;i++)
        {
         for(int j=1;j<=nlocaly;j++)
            {
              AP3=(1.0/6.0)*(P2(i+1,j)+P2(i-1,j) + P2(i,j+1) + P2(i,j-1));
              BP3=(1.0/12.0)*(P2(i-1,j-1) + P2(i-1,j+1) + P2(i+1,j-1) + P2(i+1,j+1));
              CP4(i,j)=AP3 + BP3;
	    //  printf("CP4(%d,%d)=%lf\n",i,j,CP4(i,j));
           }
 
       }


 }

double Polymer2D::g(double phi)
{
  double q(0.0);
  q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
  return q;
}


void Polymer2D::SetSecondLaplacian2(MPI_Comm new_comm)
{
   MPI_Comm_rank(new_comm, &my2drank);
   MPI_Comm_size(new_comm, &size);
   // ExchangeData(new_comm,gamma_p);
   double AP=0.0, BP=0.0;
   for(int i=1;i<=nlocalx;i++)
      {   
       for(int j=1;j<=nlocaly;j++)
          {   
           AP=(1.0/6.0)*(gamma(i+1,j)    + gamma(i-1,j)
	               + gamma(i,j+1)    + gamma(i,j-1));
	   BP=(1.0/12.0)*(gamma(i-1,j+1) + gamma(i-1,j-1)
                        +gamma(i+1,j+1)  + gamma(i+1,j-1));
           Laplacian2(i,j) = AP + BP;
          }
      }

}
void Polymer2D::UpdateSolution(MPI_Comm new_comm)
{
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);      
    for(int i=1;i<nlocalx+1;i++)
       {
        for(int j=1;j<nlocaly+1;j++)
           {
            PHI_old[i][j]=PHI[i][j];
	    P3[i][j]=PHI_old[i][j];
           }
       }

}

void Polymer2D::FiniteDifferenceScheme(MPI_Comm new_comm)
{
  MPI_Comm_rank(new_comm, &my2drank);
  MPI_Comm_size(new_comm, &size);
    //double r=0.5;
  std::cout<<"in finite difference scheme parallel"<<std::endl;
  for(int i=1;i<nlocalx+1;i++)
     {
      for(int j=1;j<nlocaly+1;j++)
         {
            // PHI[i][j]= PHI_old(i,j) - B*(PHI_old(i,j)-1.0 + 2*r)+ gamma(i,j)-Laplacian2(i,j);
//	   printf("P2=[%d][%d]=%lf\n",i,j,P2[i][j]);	 
           PHI[i][j]= PHI_old(i,j)+EMME*delta_t*(-B*(1.0-PPP[i][j])*PHI_old(i,j)+ dxx*(gamma(i,j)-Laplacian2(i,j)) + dxx*CP4[i][j]);
	     PHI_local_result[(i-1)*nlocaly + j-1]= PHI[i][j];
	     //printf("CP4[%d][%d]=%lf\n",i,j,CP4[i][j]);
          }
     }

}

void Polymer2D::ExchangeData(MPI_Comm new_comm, double **array)
{
  MPI_Status status;
  MPI_Comm_rank(new_comm,&my2drank);
  int coords[2], upper_nbr(0), down_nbr(0), tag, left_nbr(0), right_nbr(0) ;
  tag=201;
  int dims[2];
  int periods[2];
  MPI_Cart_coords(new_comm,my2drank,2,coords);
  MPI_Cart_get(new_comm,2,dims,periods, coords );
  //copy boundary points for exchange among mpi processes
  if((dims[0]>1)&&(dims[1]>1))
   {
     if(dims[0]>1)
       {
         for(int j=1;j<=nlocaly; j++)
            {
             matrix_upper[j-1]=array[1][j];
            }
         for(int j=1;j<=nlocaly; j++)
	    {
              matrix_lower[j-1]=array[nlocalx][j];

            }
         MPI_Comm_rank(new_comm,&my2drank);
         MPI_Cart_coords(new_comm,my2drank,2,coords);
         MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr ); //move along x direction
         MPI_Sendrecv_replace(matrix_upper, nlocaly,  MPI_DOUBLE, upper_nbr, tag,
	          	      down_nbr,tag, new_comm,&status);
         MPI_Sendrecv_replace(matrix_lower, nlocaly,  MPI_DOUBLE, down_nbr, tag,
		           upper_nbr,tag, new_comm,&status);

      //copy received ghost points into matrix, array

         for(int j=1; j<=nlocaly; j++)
            {
                array[0][j]=matrix_lower[j-1];
		array[nlocalx+1][j]=matrix_upper[j-1];
	     }
        }
      if(dims[1]>0)
        {
         for(int i=0;i<=nlocalx+1;i++)
	    {
             matrix_left[i]=array[i][1];
	    }
         for(int i=0;i<=nlocalx+1;i++)
	   {
             matrix_right[i]=array[i][nlocaly];
	   }
         MPI_Comm_rank(new_comm,&my2drank);
         MPI_Cart_coords(new_comm,my2drank,2,coords);
         MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
         MPI_Sendrecv_replace(matrix_right, nlocaly,  MPI_DOUBLE, right_nbr, tag,
		          left_nbr,tag, new_comm,&status);
         MPI_Sendrecv_replace(matrix_left, nlocaly,  MPI_DOUBLE, left_nbr, tag,
		          right_nbr,tag, new_comm,&status);
     //copy received ghost points into matrix

        for(int i=0;i<=nlocalx+1; i++)
	   {
             array[i][0]=matrix_right[i];
	     array[i][nlocaly+1]=matrix_left[i];
	   }


        }
   }
else if((dims[0]>1)&&(dims[1]<1))
  {
     for(int j=1;j<=nlocaly; j++)
            {
             matrix_upper[j-1]=array[1][j];
            }
         for(int j=1;j<=nlocaly; j++)
	    {
              matrix_lower[j-1]=array[nlocalx][j];

            }
         MPI_Comm_rank(new_comm,&my2drank);
         MPI_Cart_coords(new_comm,my2drank,2,coords);
         MPI_Cart_shift( new_comm, 0, 1, &upper_nbr, &down_nbr ); //move along x direction
         MPI_Sendrecv_replace(matrix_upper, nlocaly,  MPI_DOUBLE, upper_nbr, tag,
	          	      down_nbr,tag, new_comm,&status);
         MPI_Sendrecv_replace(matrix_lower, nlocaly,  MPI_DOUBLE, down_nbr, tag,
		           upper_nbr,tag, new_comm,&status);

      //copy received ghost points into matrix, array

         for(int j=1; j<=nlocaly; j++)
            {
                array[0][j]=matrix_lower[j-1];
		array[nlocalx+1][j]=matrix_upper[j-1];
	     }
  }
else if((dims[0]<1)&&(dims[1]>1))
{
  std:: cout<<"Not yet implemented. Realigned the processors. That is, let Procy=1 and Procx=size"<<std::endl;
}

}


void Polymer2D::WriteToFile_MPI(MPI_Comm new_comm)
{
  MPI_Status status;
  MPI_File     fh;
  MPI_Datatype filetype;
  int dims[2],coords[2], periods[2], start_indices[2], localsizes[2];
  int globalsizes[2];

  globalsizes[0]=BASE::Nx;
  globalsizes[1]=BASE::Ny;
  localsizes[0]=(BASE::Nx/BASE::Procx); localsizes[1]=(BASE::Ny/BASE::Procy) ;
  int local_array_size =BASE::nlocalx*BASE::nlocaly;

  MPI_Cart_get (new_comm,3,dims,periods, coords );
  MPI_Comm_rank(new_comm, &my2drank);
  MPI_Comm_size(new_comm, &size);
  MPI_Cart_coords(new_comm,my2drank,2,coords);

  start_indices[0]=coords[0]*localsizes[0];
  start_indices[1]=coords[1]*localsizes[1];
  printf("rank=%d,r=%d,start_indices=(%d,%d)\n",my2drank,Nx%Procx,start_indices[0],start_indices[1]);
  MPI_Type_create_subarray(2, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  char outputfilename[]="datafile_mpi";
  char filememview[]="native";
  MPI_File_open(new_comm, outputfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, filememview,MPI_INFO_NULL);
  MPI_File_write_all(fh, PHI_local_result,local_array_size , MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  MPI_Type_free(&filetype);

}
