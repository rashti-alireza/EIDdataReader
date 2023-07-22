/* Alireza Rashti Feb 2020 */

#include "bam.h"
#include "elliptica_id_reader_lib.h"
#include "bam_EIDdataReader.h"

#define STR_LEN_MAX (1000)

/* get binary black hole neutron star initial data into bam */
void bam_EIDdataReader(void) 
{
  if(!Getv("physics", "EIDdataReader")) return;
  
  /* functions */
  AddFun(INITIALDATA_SET, EIDdataReader,
      "Import initial data made by Elliptica");
  
  /* parameters */
  printf("Adding EID initial data reader: read data from Elliptica\n");

  AddPar("EIDdataReader_exe", "./elliptica", 
         "location of elliptica executable");
  AddPar("EIDdataReader_checkpoint", "",
         "/full/path/to/elliptica/checkpoint/file");
  
  AddPar("EIDdataReader_physics", "BHNS",
  "The physical system that the initial data were constructed."
  "options:"
  "BHNS: BH-NS binary"
  "NSNS: NS-NS binary"
  "SBH : single BH");
  
  AddPar("EIDdataReader_BHfiller","ChebTn_Ylm_perfect_s2",
   "options:"
   "ChebTn_Ylm_perfect_s2:"
   "  fill excised BH with data demanding C2 continuity for perfect S2"
   "ChebTn_general_s2:"
   "  fill excised BH with data demanding C2 continuity for general S2"
   "none: no filling.");
  
  AddPar("EIDdataReader_inf_shift","zero",
   "behavior of shift at infinity."
   "options: [zero,inspiral]"
   "o. zero     = for inertial frame."
   "o. inspiral = for corotating frame.");
  
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


