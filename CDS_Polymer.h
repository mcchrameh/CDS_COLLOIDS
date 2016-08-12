#ifndef _POLYMER_
#define _POLYMER_
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "CDS_BASE.h"
#include "mpi.h"


class Polymer2D: virtual public  BASE
{
  protected:
	    int Nx, Ny;
	  //  double delta_t, M, C1,C2, M, tau, A, F, V,dxx, B, r;
	  //  double Max_time_iteration;


	  //  std::vector<std::vector<double> > PHI_p; //moved to base class
	  //  std::vector<std::vector<double> > PHI_old_p; //moved to base class
	 //   std::vector<std::vector<double> > Laplacian2_p;
	 //   std::vector<std::vector<double> > gamma_p;
	    std:: vector<int> up1x;
	    std:: vector<int> down1x;
    

	   // void Read_input_parameters(int *integers, double *rationals);



 public:
	Polymer2D();
        ~Polymer2D();
        void UpdateSolution(MPI_Comm new_comm);
        // void initialCondition();
        void setLaplacianBase(MPI_Comm new_comm);
        double g(double phi);
        void SetSecondLaplacian2(MPI_Comm new_comm);
        void FiniteDifferenceScheme(MPI_Comm new_comm);
        void ExchangeData(MPI_Comm new_comm, double **array);
	
    //	void Solver();
	std::string make_output_filename(int index)
	    {
	      std::ostringstream ss;
              ss << "result_" << index << ".vtk";
	      return ss.str();
	    }
	std::string make_output_filename1(int index)
	   {
	      std::ostringstream ss1;
              ss1 << "output1_" << index << ".dat";
              return ss1.str();
           }

         void WriteToFile_MPI(MPI_Comm new_comm);
	 void ExchangeData(MPI_Comm new_comm, std::vector<std::vector<double> > **array);
	 MPI_Comm CreatCartesianTopology();
//	 void FreeCommunicator(MPI_Comm new_comm);










};
#endif
