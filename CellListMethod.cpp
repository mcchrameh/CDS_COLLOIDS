#include <isostream>








void GenerateColloids()
{
  double c[3],gap[3],e[3],vSum[3],gvSum[3],vMag;
  int j,a,nX,nY,nZ;
  double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
	                   {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}; 

  
 /* Set up a face-centered cubic (fcc) lattice */
for(a=0; a<3; a++) gap[a] = al[a]/InitUcell[a];
    n = 0;
   for(nZ=0; nZ<InitUcell[2]; nZ++) 
      {
       c[2] = nZ*gap[2];
       for(nY=0; nY<InitUcell[1]; nY++) 
          {
	   c[1] = nY*gap[1];
	   for(nX=0; nX<InitUcell[0]; nX++) 
	      {
	       c[0] = nX*gap[0];
	       for (j=0; j<4; j++) 
	           {
		    for(a=0; a<3; a++)
		       r[n][a] = c[a] + gap[a]*origAtom[j][a];
										 
		       ++n;
		    }
	       }
	   }
      }
  /* Total # of atoms summed over processors */
  MPI_Allreduce(&n,&nglob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if(sid == 0) printf("nglob = %d\n",nglob);


}
