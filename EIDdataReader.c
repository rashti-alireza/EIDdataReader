/* Alireza Rashti Feb 2020 */

#include "bam.h"
#include "elliptica_id_reader_lib.h"
#include "EIDdataReader.h"

/* initialize the following fields imported using the initial data.
// NOTE: in Elliptica fields are the same as variables in BAM.
// NOTE: for more info about the following fields look:
// 'bam_adm.c' and 'bam_grhd.c' files.
// NOTE: we assume the upper and lower inidices defined in those above
// files won't be changing, otherwise we need to change them in Elliptica 
// as well.
// NOTE: initial data generally are made in a corotating frame. 
// when importing the initial data in BAM, there must be a transformation
// from this corotating frame to another frame such that the new one 
// asymptotically behaves like an inertial frame. In practice, 
// we ignore omega *r part of shifts vector and transform other fields
// consequently.  */


/* algorithm:
// ==========
// at each level:
// 1. write grid points in Cartesian coords. (x,y,z) into a file
// 2. call Elliptica and read this file composed of (x,y,z) coords.
// 3. interpolate the required fields by bam at (x,y,z)
// 4. write the values of interpolated fields into another file.
// 5. read this file by bam to initializing the fields. 
// ->return value: 0 */
int EIDdataReader(tL *const level)
{
  printf("{ Importing initial data from Elliptica ...\n");
  
  const char *const outdir = Gets("outdir");
  //const int rank = bampi_rank();
  //char coords_file_path[STR_LEN_MAX] = {'\0'};
  //char fields_file_path[STR_LEN_MAX] = {'\0'};
  char checkpt_path[STR_LEN_MAX];
  sprintf(checkpt_path, "%s/%s", Gets("EIDdataReader_datadir"), "checkpoint.dat");
  const char *checkpnt_path  = checkpt_path;
  int i;

  /* read fields from the file made by Elliptica */
  FILE *file = 0;
  file = fopen(checkpnt_path,"r");
  if(!file) 
    errorexits("could not open %s", checkpnt_path);

  //const double *const x = level->v[Ind("x")];
  //const double *const y = level->v[Ind("y")];
  //const double *const z = level->v[Ind("z")];

  /* initialize ID Reader */ 
  Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path,"generic");
  idr->npoints  = level->npoints;
  idr->x_coords = level->v[Ind("x")];
  idr->y_coords = level->v[Ind("y")];
  idr->z_coords = level->v[Ind("z")];

  /* read fields content */
  if(Getv("EIDdataReader_physics",BHNS))
  {
    idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,grhd_rho,grhd_p,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
    idr->set_param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
    idr->set_param("ADM_B1I_form","zero",idr); 
  }
  else if(Getv("EIDdataReader_physics",SBH))
  {
    idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz";
  }
  else
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  
  /* interpolate */
  elliptica_id_reader_interpolate(idr);
  
  printf("~> Populating BAM variables based on initial data ...\n");
  fflush(stdout);
  forallpoints(level, i)
    {
      Ptr(level, "alpha")[i] = idr->field[idr->indx("alpha")][i];
      Ptr(level, "betax")[i] = idr->field[idr->indx("betax")][i];
      Ptr(level, "betay")[i] = idr->field[idr->indx("betay")][i];
      Ptr(level, "betaz")[i] = idr->field[idr->indx("betaz")][i];

      Ptr(level, "adm_gxx")[i] = idr->field[idr->indx("adm_gxx")][i];
      Ptr(level, "adm_gxy")[i] = idr->field[idr->indx("adm_gxx")][i];
      Ptr(level, "adm_gxz")[i] = idr->field[idr->indx("adm_gxx")][i];
      Ptr(level, "adm_gyy")[i] = idr->field[idr->indx("adm_gxx")][i];
      Ptr(level, "adm_gyz")[i] = idr->field[idr->indx("adm_gxx")][i];
      Ptr(level, "adm_gzz")[i] = idr->field[idr->indx("adm_gxx")][i];

      Ptr(level, "adm_Kxx")[i] = idr->field[idr->indx("adm_Kxx")][i];
      Ptr(level, "adm_Kxy")[i] = idr->field[idr->indx("adm_Kxx")][i];
      Ptr(level, "adm_Kxz")[i] = idr->field[idr->indx("adm_Kxx")][i];
      Ptr(level, "adm_Kyy")[i] = idr->field[idr->indx("adm_Kxx")][i];
      Ptr(level, "adm_Kyz")[i] = idr->field[idr->indx("adm_Kxx")][i];
      Ptr(level, "adm_Kzz")[i] = idr->field[idr->indx("adm_Kxx")][i];

      if(Getv("EIDdataReader_physics",BHNS)){
        Ptr(level, "grhd_rho")[i] = idr->field[idr->indx("grhd_rho")][i];
        Ptr(level, "grhd_epsl")[i] = idr->field[idr->indx("grhd_epsl")][i];
        Ptr(level, "grhd_p")[i] = idr->field[idr->indx("grhd_p")][i];
        Ptr(level, "grhd_vx")[i] = idr->field[idr->indx("grhd_vx")][i];
        Ptr(level, "grhd_vy")[i] = idr->field[idr->indx("grhd_vy")][i];
        Ptr(level, "grhd_vz")[i] = idr->field[idr->indx("grhd_vz")][i];
      }
    }

  
  printf("~> Populating BAM variables based on initial data --> Done.\n");
  fflush(stdout);
  
  elliptica_id_reader_free(idr);
  printf("} Importing initial data from Elliptica --> Done.\n");
  fflush(stdout);
  return 0;
}
