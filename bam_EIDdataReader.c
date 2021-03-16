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
  
  /* set output directory */
  AddPar("EIDdataReader_outdir", "elliptica_id",
         "location of imported initial data for bam");

  AddPar("EIDdateReader_physics", "BHNS",
  "The physical system which this initial data was constructed."
  "options:"
  "BHNS: BH-NS binary"
  "SBH : single BH");
  
  AddPar("EIDdateReader_BHfiller","ChebTn_Ylm_perfect_s2",
   "options:"
   "ChebTn_Ylm_perfect_s2:"
   "  fill excised BH with data demanding C2 continuity for perfect S2"
   "ChebTn_general_s2:"
   "  fill excised BH with data demanding C2 continuity for general S2"
   "none: no filling.");
  
  AddPar("EIDdateReader_x_CM","0","elliptica's x_CM of system");
  AddPar("EIDdateReader_y_CM","0","elliptica's y_CM of system");
  AddPar("EIDdateReader_z_CM","0","elliptica's z_CM of system");
  AddPar("mass1","1.0", "used in grid setup at PRE_GRID");
  AddPar("mass2","1.0", "used in grid setup at PRE_GRID");
  
}


