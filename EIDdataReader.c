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
  const int rank = bampi_rank();
  char coords_file_path[STR_LEN_MAX] = {'\0'};
  char fields_file_path[STR_LEN_MAX] = {'\0'};
  
  /* if this is a checkpoint, just read field values from files */
  if(Getv("checkpoint","restart"))
  {
    sprintf(fields_file_path, "%s_previous/%s/fields_level%d_proc%d.dat",
      outdir,Gets("EIDdataReader_outdir"), level->l, rank);

    /* populate fields for bam */
    populate_fields_for_bam(level,fields_file_path);
  }
  /* if requested to load from specific directory */
  else if (!Getv("EIDdataReader_loadfrom","not_set!"))
  {
    sprintf(fields_file_path, "%s/fields_level%d_proc%d.dat",
        Gets("EIDdataReader_loadfrom"), level->l, rank);
    
    /* populate fields for bam */
    populate_fields_for_bam(level,fields_file_path);
  }
  else/* make everything from scratch */
  {    
    /* populate fields for bam */
    populate_fields_for_bam(level,Gets("EIDdataReader_datadir"));
  }
  
  printf("} Importing initial data from Elliptica --> Done.\n");
  fflush(stdout);
  return 0;
}

/* populating the fields made by Elliptica into bam */
static void populate_fields_for_bam(tL *const level, char *const checkpoint_dir)
{
  printf("~> Populating BAM variables based on initial data ...\n");
  fflush(stdout);
  FILE *file = 0;
  double *v  = 0;
  char *v_name = 0;
  char checkpt_path[STR_LEN_MAX];
  sprintf(checkpt_path, "%s/%s", checkpoint_dir, "checkpoint.dat");
  const char *checkpnt_path  = checkpt_path;
  const double *const x = level->v[Ind("x")];
  const double *const y = level->v[Ind("y")];
  const double *const z = level->v[Ind("z")];
  
  char *match_str;
  char msg[STR_LEN_MAX];
  int msg_len = (int)strlen(END_MSG)+1;
  int i,f;
  
  /* read fields from the file made by Elliptica */
  file = fopen(checkpnt_path,"r");
  if(!file) 
    errorexits("could not open %s", checkpnt_path);
    
  // initialize ID Reader
  Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path,"generic");
  idr->npoints  = level->npoints;
  idr->x_coords = x;
  idr->y_coords = y;
  idr->z_coords = z;

  /* read fields content */
  if(Getv("EIDdataReader_physics",BHNS))
  {
    idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,grhd_rho,grhd_p,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
  }
  else if(Getv("EIDdataReader_physics",SBH))
  {
    idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz";
  }
  else
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  
  idr->set_param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
  idr->set_param("ADM_B1I_form","zero",idr);  

  // interpolate
  elliptica_id_reader_interpolate(idr);

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

  /*double *gxx = idr->field[idr->indx("adm_gxx")];//Ptr(level, "adm_gxx");
  double *gxy = idr->field[idr->indx("adm_gxy")];//Ptr(level, "adm_gxy");
  double *gxz = idr->field[idr->indx("adm_gxz")];//Ptr(level, "adm_gxz");
  double *gyy = idr->field[idr->indx("adm_gyy")];//Ptr(level, "adm_gyy");
  double *gyz = idr->field[idr->indx("adm_gyz")];//Ptr(level, "adm_gyz");
  double *gzz = idr->field[idr->indx("adm_gzz")];//Ptr(level, "adm_gzz");*/

  
  /* set tr(K0)? */
 /* if (0)
  {
    i = Ind("adm_K0"); enablevar(level, i);
    double *Kxx = idr->field[idr->indx("adm_Kxx")];//Ptr(level, "adm_Kxx");
    double *Kxy = idr->field[idr->indx("adm_Kxy")];//Ptr(level, "adm_Kxy");
    double *Kxz = idr->field[idr->indx("adm_Kxz")];//Ptr(level, "adm_Kxz");
    double *Kyy = idr->field[idr->indx("adm_Kyy")];//Ptr(level, "adm_Kyy");
    double *Kyz = idr->field[idr->indx("adm_Kyz")];//Ptr(level, "adm_Kyz");
    double *Kzz = idr->field[idr->indx("adm_Kzz")];//Ptr(level, "adm_Kzz");
    double *K0  = idr->field[idr->indx("adm_K0")];//Ptr(level, "adm_K0");
    
    double ixx, ixy, ixz, iyy, iyz, izz;
    
    forallpoints(level, i)
    {
      invg(gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
           &ixx, &ixy, &ixz, &iyy, &iyz, &izz);

      K0[i] = 
          (ixx*Kxx[i] + iyy*Kyy[i] + izz*Kzz[i] 
              + 2*(ixy*Kxy[i] + ixz*Kxz[i] + iyz*Kyz[i]));

    }
  }*/
  
  /* some tests */
  //double detg;
  
  /* check if determinant is positvie, since the metric is Reimannian */
  /*forallpoints(level,i)
  {
    detg=(2.*gxy[i]*gxz[i]*gyz[i] + gxx[i]*gyy[i]*gzz[i] -
             gzz[i]*gxy[i]*gxy[i] - gyy[i]*gxz[i]*gxz[i] -
             gxx[i]*gyz[i]*gyz[i]);
    
    if(detg<=0)
    {
      printf("det(g_ij)=%g at ccc=i=%d:  x=%g y=%g z=%g\n", 
             detg, i, x[i], y[i], z[i]);
    }
  }*/

  elliptica_id_reader_free(idr);
  printf("~> Populating BAM variables based on initial data --> Done.\n");
  fflush(stdout);
}
