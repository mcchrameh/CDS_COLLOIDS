#ifndef _COLLIOD_
#define _COLLIOD_
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "CDS_BASE.h"

#include "mpi.h"

/*
typedef struct
{
 int body_index;	
double x;
double y;
double vx;
double vy;

}Body;
*/
class Colloid:  virtual public BASE
{
  protected:

	  // double ENNE, ALPHA, GAMMs, GAMMb ;
	    std::vector<double> Posx;
            std::vector<double> Posy;
	    std::vector<double> Vx;
	    std::vector<double> Vy;
	   
	    std::vector<double> Fx;
	    std::vector<double> Fy;
           
	   
            void ComputeForce(MPI_Comm new_comm);//, BODY* bd, Node* root, double dist);
	   
	    void UpdatePosition();


  public:
        Colloid();
        ~Colloid();
       void Exchange_BcColloids(MPI_Comm new_comm);
      // void WriteToFile_Colloids();



};
#endif
