#include "CDS_Colloids.h"
#include <random>
#include <time.h>


#define Random_min  -0.05
#define Random_max  0.05 








Colloid::Colloid():BASE()
        {
          Posx.resize(Nparticles);
	  Posy.resize(Nparticles);
	  Fx.resize(Nparticles);
	  Fy.resize(Nparticles);
	  Vx.resize(Nparticles);
	  Vy.resize(Nparticles);



	}
Colloid::~Colloid()
{

 printf("Destructor of Colloid class called\n");

}
void Colloid::Exchange_BcColloids(MPI_Comm new_comm)
{
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);



 





}
void Colloid::ComputeForce(MPI_Comm new_comm)//, BODY* bd, Node* root, double dist=0.0)
{
 MPI_Comm_rank(new_comm, &my2drank);
 MPI_Comm_size(new_comm, &size);
  double G4(0.0), AA(0.0), BB(0.0), NEW1(0.0), NEW2(0.0), G2(0.0), G3(0.0), G3A(0.0),C(0.0);
  double DSUMX(0.0),DSUMY(0.0), DSUMXY(0.0);
/*---Make a linked cell list, LinkedCells-----*/
  int Ncells;
  Ncells = LinkedCells[0]*LinkedCells[1];//<--total number of cells on a process grid
  std::vector<int> mc;//<---use for calculating the vector cell index to which particle i belongs
  mc.resize(2);
  std::vector<int> mc1;
  mc1.resize(2);
  int c(0), c1(0),i,j;
  for(c=0;c<Ncells;c++)
     {
	head[c]=EMPTY;
     }
  /* Scan atoms to construct headers, head, & linked lists, lscl */
    
    for ( auto i=0; i<(int)bd.size(); i++)
        {
	  for (int direction=0; direction<2; direction++)
	      {
		//printf("bd[%d].r[%d]=%lf\n",i, direction,bd[i].r[direction]);
		
		mc.at(direction) = (int)((bd[i].r[direction] + CellSize[direction])/CellSize[direction]);
		//printf("mc[%d]=%lf\n",direction,CellSize[direction]);
	      }
          /* Translate the vector cell index, mc, to a scalar cell index */
	   c = mc[0]*LinkedCells[1] +mc[1];
          /* Link to the previous occupant (or EMPTY if you're the 1st) */
           LinkedCell_list[i] = head[c];
          /* The last one goes to the header */
	   head[c] = i;
        } /* Endfor atom i */

    // Scan inner cells /
      for (mc[0]=1; mc[0]<=LinkedCells[0]; (mc[0])++)
	  {
            for (mc[1]=1; mc[1]<=LinkedCells[1]; (mc[1])++)
	         {
	           /// Calculate a scalar cell index /
		     c = mc[0]*LinkedCells[1] + mc[1];
		   // Skip this cell if empty /
		    if (head[c] == EMPTY) continue;
		    // Scan the neighbor cells (including itself) of cell c /
		     for(mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
			  {
		           for(mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++)
		                {
                                   // Calculate the scalar cell index of the neighbor cell /
				   c1 = mc1[0]*LinkedCells[1] + mc1[1];
				    // Skip this neighbor cell if empty /
				    if (head[c1] == EMPTY) continue;
                                     // Scan atom i in cell c /
				     i = head[c];
				     while (i != EMPTY) 
				     {
					// Scan atom j in cell c1 /
					j = head[c1];
				        while (j != EMPTY) 
					{
					 // No calculation with itself /
				          if (j != i) 
					    {
                                              double dx=bd[i].r[0]-bd[j].r[0];
                                              double dy=bd[i].r[1]-bd[j].r[1];
                                             //dx=dx-nlocalx*(int((dx/(nlocalx/2)+3)/2)-1.0);
                                            //dy=dy-nlocaly*(int((dy/(nlocaly/2)+3)/2)-1.0);
                                             double R2=sqrt(pow(dx,2)+pow(dy,2));
                                            if(R2<R12[i][j])
                                              {
                                               Dx[i][j]=dx;
	                                       Dy[i][j]=dy;
                                               RI[i][j]=R2;
	                                       RI[j][i]=R2;
                                               C = U1/R12[i][j];
                                               AA=RI[i][j]/R12[i][j];
	                                       G2 = pow(AA + beta,2.0);
	                                       G3 = 1.0 + ALPHA*(AA + beta);
	                                       BB = (AA-1.0)*(-ALPHA);
	                                       G3A = exp(BB);
	                                       NEW1=pow((1.0 + (ENNE/ALPHA)+ beta),2.0);
	                                       NEW2=(exp(-ENNE))*((1.0 + ALPHA*(beta + (ENNE/ALPHA)+1.0)))/(NEW1);
	                                       G4 = ((C*(G3A)*G3)/G2)-(C*NEW2);

                                               Fx[i]=Fx[i] + G4*(Dx[i][j]/RI[i][j]);
	                                       Fy[i]=Fy[i] + G4*(Dy[i][j]/RI[i][j]);
	                                       //   printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);
                                               // Fx[j]=Fx[j] + G4*(Dx[j][i])/RI[i][j];
                                               // Fy[j]=Fy[j] + G4*(Dy[j][i])/RI[i][j];
                                               }
                                           }
                                       
		                           j = LinkedCell_list[j];
                                         }
				         i = LinkedCell_list[i];
                                        }
			      }
			  }
			}
		  }
  //Compute random force
  double csi1(0.0), csi2(0.0);
  srand(time(NULL) + rank);
  double GAMM(0.0);




  for(int i=N_start;i<N_start + nlocal_particles; i++)
 // for(int i=0;i< nlocal_particles; i++)
     {
       double range =Random_max-Random_min;
       double div =RAND_MAX / range;
       csi1 = Random_min + (rand()/div);
       csi2 = Random_min + (rand()/div);
      //  printf("Random number (csi1=%lf, csi2=%lf) \n", csi1, csi2);

       if(bd[i].index==0)
        {
	  GAMM= GAMMs;
          bd[i].v[0]= (1.0/GAMM)*(Fx[i]+ TOTX[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
          bd[i].v[1]= (1.0/GAMM)*(Fy[i]+ TOTY[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi1*sqrt(delta_t);
	 // printf("Colloid Force (fx=%lf, fy=%lf) \n", Fx[i], Fy[i]);

        }
       if(bd[i].index==1)
	{
	  GAMM= GAMMb;
          bd[i].v[0]= (1.0/GAMM)*(Fx[i]+ TOTX[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
          bd[i].v[1]= (1.0/GAMM)*(Fy[i]+ TOTY[i])*delta_t + sqrt(2.0*temp1/GAMM)*csi2*sqrt(delta_t);
         // printf("Colloid (TOTX =%lf, TOTY=%lf) \n", TOTX[i], TOTY[i]);
	}
     //position update
     bd[i].r[0] =bd[i].r[0] + bd[i].v[0]; //<---bd[i].vx is already multiplied by t. 
     bd[i].r[1] =bd[i].r[1] + bd[i].v[1];
   //  printf("bd[%d].x=%lf, bd[%d].y=%lf\n",i,bd[i].x,i,bd[i].y);
      //bc condition
     //bd[i].x=bd[i].x-nlocalx*((int)(bd[i].x/nlocalx + 1.0)-1.0);
     if((bd[i].r[0]>(double)nlocalx)||(bd[i].r[1]>(double)nlocaly))
	{
          //std::vector<BODY>::iterator it;
          printf("bd[%d].x=%lf, bd[%d].y=%lf\n",i,bd[i].r[0],i,bd[i].r[1]);		
	  send_bodies.push_back(bd[i]);
//	  bd.erase(bd.begin()+i-1);

	}
      //bd[i].x=bd[i].x-(nlocalx/2)*((double)rand())/(RAND_MAX);


     //printf("Rank=%d,Positionx=%lf\n ", my2drank, bd[i].x);
     //bd[i].y=bd[i].y-nlocaly*((int)(bd[i].y/nlocaly + 1.0)-1.0);

     //computing mean displacement
     Posx[i]=Posx[i] + bd[i].v[0];
     Posy[i]=Posy[i] + bd[i].v[1];

     DSUMX  = DSUMX + pow(Posx[i],2.0);
     DSUMY  = DSUMY + pow(Posy[i],2.0);
     DSUMXY = DSUMXY + Posx[i]*Posy[i];
    
    // printf("Colloid positions (x=%lf, y=%lf) \n", bd[i].x, bd[i].y);



    }
 // DSUMX=DSUMX/

  
}

/*
void Colloid::WriteToFile_colloid(MPI_Comm new_comm)
{
   MPI_Status status;
   MPI_File     fh;
   MPI_Datatype filetype;
   int dims[2],coords[2], periods[2], start_indices[2], localsizes[2];
   int globalsizes[2];
   

}*/

