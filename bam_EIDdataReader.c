/* Alireza Rashti Feb 2020 */

#include "bam.h"
#include "bam_EIDdataReader.h"

/* get binary black hole neutron star initial data into bam */
void bam_EIDdataReader(void) 
{
  if(!Getv("physics", "EIDdataReader")) return;
  
  /* functions */
  AddFun(PRE_GRID, EIDpreGrid,
      "set some pars important for making the grid");
  AddFun(INITIALDATA_SET, EIDdataReader,
      "Import initial data made by Elliptica");
  
  /* parameters */
  printf("Adding EID initial data reader: read data from Elliptica\n");

  AddPar("EIDdataReader_exe", "./elliptica", 
         "location of elliptica executable");
  AddPar("EIDdataReader_datadir", "",
         "location of elliptica outdir with data");
  
  AddPar("EIDdateReader_physics", "BHNS",
  "The physical system which this initial data was constructed."
  "options:"
  "BHNS: BH-NS binary"
  "SBH : single BH");
  
  AddPar("EIDdateReader_BHfiller","ChebTn_Ylm",
   "ChebTn_Ylm:fill excised BH with data demanding C2 continuity "
   "across horizon and extrapolant is "
   "f(r,th,ph) = C_{ilm}*ChebyshevT(i,r)*Ylm(th,ph).");
  
  AddPar("EIDdateReader_x_CM","0","elliptica's x_CM of system");
  AddPar("EIDdateReader_y_CM","0","elliptica's y_CM of system");
  AddPar("EIDdateReader_z_CM","0","elliptica's z_CM of system");
  AddPar("mass1","1.0", "used in grid setup at PRE_GRID");
  AddPar("mass2","1.0", "used in grid setup at PRE_GRID");
  
}


