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
static const char *const import_fields_with_matter[] = /* matter included */
{
"alpha",/* lapse: alpha */
"betax","betay","betaz",/* shift: beta^i */
"adm_gxx","adm_gxy","adm_gxz",/* metric: g_ij */
"adm_gyy","adm_gyz","adm_gzz",/* metric: g_ij */
"adm_Kxx","adm_Kxy","adm_Kxz",/* extrinsic curvature: K_ij */
"adm_Kyy","adm_Kyz","adm_Kzz",/* extrinsic curvature: K_ij */

/* matter part */
"grhd_rho",/* primitive rho */
"grhd_p",/* primitive p */
"grhd_epsl",/* primitive epsilon: total_energy_density = grhd_rho(1+grhd_epsl)*/
"grhd_vx","grhd_vy","grhd_vz",/* primitive v, measured by an Eulerian observer, 
                              // v^i = u^i/(alpha u^0) + beta^i / alpha
                              // where u^{mu}=(u^0,u^i) is the 4-velocity of the fluid. */
                         
  0/* --> detemine the last pointer */
};

static const char *const import_fields_no_matter[] = /* matter excluded */
{
"alpha",/* lapse: alpha */
"betax","betay","betaz",/* shift: beta^i */
"adm_gxx","adm_gxy","adm_gxz",/* metric: g_ij */
"adm_gyy","adm_gyz","adm_gzz",/* metric: g_ij */
"adm_Kxx","adm_Kxy","adm_Kxz",/* extrinsic curvature: K_ij */
"adm_Kyy","adm_Kyz","adm_Kzz",/* extrinsic curvature: K_ij */
  0/* --> detemine the last pointer */
};

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
  //const double *const x = level->v[Ind("x")];
  //const double *const y = level->v[Ind("y")];
  //const double *const z = level->v[Ind("z")];
  
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
  idr->x_coords = level->v[Ind("x")];
  idr->y_coords = level->v[Ind("y")];
  idr->z_coords = level->v[Ind("z")];

  /* read fields content */
  if(Getv("EIDdataReader_physics",BHNS))
  {
    idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,grhd_rho,grhd_p,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";

    const double x_CM = idr->get_param_dbl("BHNS_x_CM",idr);
    const double y_CM = idr->get_param_dbl("BHNS_y_CM",idr);
    const double z_CM = idr->get_param_dbl("BHNS_z_CM",idr);
    const double ns_center_x = idr->get_param_dbl("NS_center_x",idr);
    const double ns_center_y = idr->get_param_dbl("NS_center_y",idr);
    const double ns_center_z = idr->get_param_dbl("NS_center_z",idr);
    const double bh_center_x = idr->get_param_dbl("BH_center_x",idr);
    const double bh_center_y = idr->get_param_dbl("BH_center_y",idr);
    const double bh_center_z = idr->get_param_dbl("BH_center_z",idr);
    const double ns_m = idr->get_param_dbl("NS_baryonic_mass_current",idr);
    const double bh_m = idr->get_param_dbl("BH_irreducible_mass_current",idr);

    /* set CM params for translation */
    MySetd("EIDdataReader_x_CM",x_CM);
    MySetd("EIDdataReader_y_CM",y_CM);
    MySetd("EIDdataReader_z_CM",z_CM);

    /* set pre grid params 
 *     // note: r_elliptica = r_CM + r_bam 
 *         //       = > r_bam = r_elliptica - r_CM. */
    MySetd("mass1", ns_m);
    MySetd("mass2", bh_m);
    MySetd("px1", ns_center_x-x_CM);
    MySetd("py1", ns_center_y-y_CM);
    MySetd("pz1", ns_center_z-z_CM);
    MySetd("px2", bh_center_x-x_CM);
    MySetd("py2", bh_center_y-y_CM);
    MySetd("pz2", bh_center_z-z_CM);
    
    /* for BHfiller */
    if (1)
    {
      /* note: the notation different here and object 1 is bh. */
      MySetd("bhmass1",bh_m);
      MySetd("bhmass2",0.);
      MySetd("bhx1", bh_center_x-x_CM);
      MySetd("bhy1", bh_center_y-y_CM);
      MySetd("bhz1", bh_center_z-z_CM);
    }
  }
  else if(Getv("EIDdataReader_physics",SBH))
  {
    idr->ifields  = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz";
    // TODO
  }
  else
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  
  idr->set_param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
  idr->set_param("ADM_B1I_form","zero",idr);  

  // interpolate
  elliptica_id_reader_interpolate(idr);

  double *gxx = idr->field[idr->indx("adm_gxx")];//Ptr(level, "adm_gxx");
  double *gxy = idr->field[idr->indx("adm_gxy")];//Ptr(level, "adm_gxy");
  double *gxz = idr->field[idr->indx("adm_gxz")];//Ptr(level, "adm_gxz");
  double *gyy = idr->field[idr->indx("adm_gyy")];//Ptr(level, "adm_gyy");
  double *gyz = idr->field[idr->indx("adm_gyz")];//Ptr(level, "adm_gyz");
  double *gzz = idr->field[idr->indx("adm_gzz")];//Ptr(level, "adm_gzz");

  
  /* set tr(K0)? */
  if (0)
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
  }
  
  /* some tests */
  double detg;
  
  /* check if determinant is positvie, since the metric is Reimannian */
  forallpoints(level,i)
  {
    detg=(2.*gxy[i]*gxz[i]*gyz[i] + gxx[i]*gyy[i]*gzz[i] -
             gzz[i]*gxy[i]*gxy[i] - gyy[i]*gxz[i]*gxz[i] -
             gxx[i]*gyz[i]*gyz[i]);
    
    if(detg<=0)
    {
      printf("det(g_ij)=%g at ccc=i=%d:  x=%g y=%g z=%g\n", 
             detg, i, idr->x_coords[i], idr->y_coords[i], idr->z_coords[i]);
    }
  }

  elliptica_id_reader_free(idr);
  printf("~> Populating BAM variables based on initial data --> Done.\n");
  fflush(stdout);
}

/* set some pars important for making the grid */
int EIDpreGrid(tL *const level)
{
  /* set pre grid parameters */
  printf("Setting some parameters relevant for setting up bam's grid:\n");
  
  if(Getv("EIDdataReader_physics",BHNS))
  {
    FILE *file = 0;
    char *const file_NAME = BHNS_"properties.txt";
    double ns_m,bh_m;/* gravitational masses of NS and BH */
    double ns_center_x,ns_center_y,ns_center_z;/* NS center coords */
    double bh_center_x,bh_center_y,bh_center_z;/* BH center coords */
    double x_CM,y_CM,z_CM;/* elliptica's CM of system */
    char id_outdir[STR_LEN_MAX] = {'\0'};
    char file_path[STR_LEN_MAX] = {'\0'};
    char str[STR_LEN_MAX] = {'\0'};
    int ret;
    
    /* open ?_properties.txt */
    sprintf(str,"%s",Gets("EIDdataReader_datadir"));/* --> a/b/c */
    if (str[strlen(str)-1] == '/')/* if there is / at the end of path */
      str[strlen(str)-1] = '\0';/* if a/b/c/ */
    sprintf(id_outdir,"%s",str);
    sprintf(file_path,"%s/%s",id_outdir,file_NAME);
    file = fopen(file_path,"r");
    if(!file) 
      errorexits("could not open %s", file_path);

    /* read parameters from ?_properties.txt */
    READ_PARAMETER_FROM_FILE(x_CM,BHNS_"x_CM")
    READ_PARAMETER_FROM_FILE(y_CM,BHNS_"y_CM")
    READ_PARAMETER_FROM_FILE(z_CM,BHNS_"z_CM")
    
    /* set CM params for translation */
    MySetd("EIDdataReader_x_CM",x_CM);
    MySetd("EIDdataReader_y_CM",y_CM);
    MySetd("EIDdataReader_z_CM",z_CM);
    
    READ_PARAMETER_FROM_FILE(ns_center_x,"NS_center_x")
    READ_PARAMETER_FROM_FILE(ns_center_y,"NS_center_y")
    READ_PARAMETER_FROM_FILE(ns_center_z,"NS_center_z")
    
    READ_PARAMETER_FROM_FILE(bh_center_x,"BH_center_x")
    READ_PARAMETER_FROM_FILE(bh_center_y,"BH_center_y")
    READ_PARAMETER_FROM_FILE(bh_center_z,"BH_center_z")
    
    READ_PARAMETER_FROM_FILE(ns_m,"NS_baryonic_mass_current")
    READ_PARAMETER_FROM_FILE(bh_m,"BH_irreducible_mass_current")
    
    /* set pre grid params 
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    MySetd("mass1", ns_m);
    MySetd("mass2", bh_m);
    MySetd("px1", ns_center_x-x_CM);
    MySetd("py1", ns_center_y-y_CM);
    MySetd("pz1", ns_center_z-z_CM);
    MySetd("px2", bh_center_x-x_CM);
    MySetd("py2", bh_center_y-y_CM);
    MySetd("pz2", bh_center_z-z_CM);
    
    /* for BHfiller */
    if (1)
    {
      /* note: the notation different here and object 1 is bh. */
      MySetd("bhmass1",bh_m);
      MySetd("bhmass2",0.);
      MySetd("bhx1", bh_center_x-x_CM);
      MySetd("bhy1", bh_center_y-y_CM);
      MySetd("bhz1", bh_center_z-z_CM);
    }
    
    fclose(file);
  }
  else if(Getv("EIDdataReader_physics",SBH))
  {
    FILE *file = 0;
    char *const file_NAME = SBH_"properties.txt";
    double ns_m,bh_m;/* gravitational masses of NS and BH */
    double ns_center_x,ns_center_y,ns_center_z;/* NS center coords */
    double bh_center_x,bh_center_y,bh_center_z;/* BH center coords */
    double x_CM,y_CM,z_CM;/* elliptica's CM of system */
    char id_outdir[STR_LEN_MAX] = {'\0'};
    char file_path[STR_LEN_MAX] = {'\0'};
    char str[STR_LEN_MAX] = {'\0'};
    int ret;
    
    /* open ?_properties.txt */
    sprintf(str,"%s",Gets("EIDdataReader_datadir"));/* --> a/b/c */
    if (str[strlen(str)-1] == '/')/* if there is / at the end of path */
      str[strlen(str)-1] = '\0';/* if a/b/c/ */
    sprintf(id_outdir,"%s",str);
    sprintf(file_path,"%s/%s",id_outdir,file_NAME);
    file = fopen(file_path,"r");
    if(!file) 
      errorexits("could not open %s", file_path);

    /* read parameters from ?_properties.txt */
    READ_PARAMETER_FROM_FILE(x_CM,SBH_"x_CM")
    READ_PARAMETER_FROM_FILE(y_CM,SBH_"y_CM")
    READ_PARAMETER_FROM_FILE(z_CM,SBH_"z_CM")
    
    /* set CM params for translation */
    MySetd("EIDdataReader_x_CM",x_CM);
    MySetd("EIDdataReader_y_CM",y_CM);
    MySetd("EIDdataReader_z_CM",z_CM);
    
    READ_PARAMETER_FROM_FILE(bh_center_x,"BH_center_x")
    READ_PARAMETER_FROM_FILE(bh_center_y,"BH_center_y")
    READ_PARAMETER_FROM_FILE(bh_center_z,"BH_center_z")
    
    READ_PARAMETER_FROM_FILE(bh_m,"BH_irreducible_mass")
    
    /* set pre grid params 
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    MySetd("mass1", bh_m);
    MySetd("mass2", 0);
    MySetd("px1", bh_center_x-x_CM);
    MySetd("py1", bh_center_y-y_CM);
    MySetd("pz1", bh_center_z-z_CM);
    
    fclose(file);
  }
  else
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  
  return 0;
}
