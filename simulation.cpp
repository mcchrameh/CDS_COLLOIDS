#include "CDS_BASE.h"
#include "CDS_Polymer.h"
#include "CDS_Colloids.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h> 
#include "mpi.h"
 

/*
 BASE  ParVer;
 Polymer2D  poly;
 Colloid Col;
 */
 class Simulate:  public Polymer2D, public Colloid
 {
   protected:
	     double t1;
	     int count1;
	     MPI_Comm new_comm;

   public:
	   Simulate():Polymer2D()//,Colloid()
		   {
                      //t=0.0; count=0;
		      new_comm = Polymer2D::BASE::CreatCartesianTopology();
		   }
           ~Simulate()
	   {
	      // Polymer2D::BASE::FreeCommunicator(new_comm);   
               std::cout<<"Destructor of simulate class called"<<std::endl;

	   }
	  void  Solver()
	  {
		  double t=0.0;
		  int count =0;

           
             
	     new_comm = Polymer2D::BASE::CreatCartesianTopology();
	    Polymer2D:: InitializePolymer(new_comm);
	    Colloid::BASE::GenerateColloids(new_comm);
           
	   //  Colloid::BASE::DistributeColloids(new_comm);
	     while(t<Max_time_iteration)
	          {
		   // std::cout<<"in while loop solver"<<std::endl;
			//  Colloid::   BASE::SetEffectiveRadius(new_comm);
			 // Colloid:: GenerateVerlet( new_comm);
			 // Colloid::BASE::DistributeColloids(new_comm);
                 	  Colloid::   BASE::SetEffectiveRadius(new_comm);
                          Polymer2D:: ExchangeData(new_comm,PHI_old);
			  Colloid::   BASE::SendBodies(new_comm);
                          Colloid::   BASE::Coupling(new_comm);
			  Colloid::   ComputeForce( new_comm);
  			  Polymer2D:: ExchangeData(new_comm,PHI_old);
			  Polymer2D:: ExchangeData(new_comm,P2);
		          Polymer2D:: setLaplacianBase( new_comm);
		          Polymer2D:: ExchangeData(new_comm,gamma);
		          Polymer2D:: SetSecondLaplacian2(new_comm);
		          Polymer2D:: ExchangeData(new_comm,Laplacian2);
		          Polymer2D:: FiniteDifferenceScheme(new_comm);
		          Polymer2D:: UpdateSolution(new_comm);
                          Colloid::   BASE::SendBodies(new_comm);

			  //Colloid::BASE::DistributeColloids(new_comm);

		     t+=1.0;
		     count++;
		     printf("time=%lf, count=%d\n",t,count);
		}
	    Polymer2D:: WriteToFile_MPI( new_comm);
            Polymer2D::BASE:: FreeCommunicator(new_comm);
           
           // MPI_Comm_free(new_comm);
         // Polymer2D::BASE::FreeCommunicator(new_comm);
	  }
	 // MPI_Comm_free(new_comm);
	 // Polymer2D::BASE::FreeCommunicator(new_comm);
 };

BASE  *Base;
Polymer2D  *Polymer2D;
Colloid *Colloid;

int main(int argc, char **argv)
{
 
 
 MPI_Init(&argc, &argv);
 Simulate code;
 code.Solver();

 MPI_Finalize();


 return 0;

}
