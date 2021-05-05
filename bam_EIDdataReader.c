/* Alireza Rashti Feb 2020 */

#include "bam.h"
#include "bam_EIDdataReader.h"

#define STR_LEN_MAX (1000)

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
  
  /* load from this dir */
  AddPar("EIDdataReader_loadfrom", "not_set!",
         "if this param is set, it won't call Elliptica anymore and "
         "just reads data from the given directory."
         "Used to quickly load data into bam for debug purposes.");
  
  /* set output directory */
  AddPar("EIDdataReader_outdir", "elliptica_id",
         "location of imported initial data for bam");

  /* mkdir outdir if not exists */
  char mydir[STR_LEN_MAX];
  sprintf(mydir,"%s/%s",Gets("outdir"),Gets("EIDdataReader_outdir"));
  if (processor0 && !system_isdir(mydir))
    if (system_mkdir(mydir)) errorexit("mkdir failed!");

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
  
  AddPar("EIDdateReader_inf_shift","zero",
   "behavior of shift at infinity."
   "options: [zero,inspiral]"
   "o. zero     = for inertial frame."
   "o. inspiral = for corotating frame.");
  
  AddPar("EIDdateReader_x_CM","0","elliptica's x_CM of system");
  AddPar("EIDdateReader_y_CM","0","elliptica's y_CM of system");
  AddPar("EIDdateReader_z_CM","0","elliptica's z_CM of system");
  AddPar("mass1","1.0", "used in grid setup at PRE_GRID");
  AddPar("mass2","1.0", "used in grid setup at PRE_GRID");
  
  AddPar("bhmass1","0.",  "used in BHfiller");
  AddPar("bhmass2","0.",  "used in BHfiller");
  AddPar("bhx1","0.",     "used in BHfiller");
  AddPar("bhy1","0.", 	  "used in BHfiller");
  AddPar("bhz1","0.", 	  "used in BHfiller");
  AddPar("bhx2","0.",     "used in BHfiller");
  AddPar("bhy2","0.",     "used in BHfiller");
  AddPar("bhz2","0.",     "used in BHfiller");
}


