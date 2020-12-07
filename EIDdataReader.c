/* Alireza Rashti Feb 2020 */

#include "bam.h"
#include "EIDdataReader.h"

/* initialize the following fields using the initial data.
// NOTE: in Elliptica fields are the same as variables in BAM.
// NOTE: for more info about the following fields look:
// 'bam_adm.c' and 'bam_grhd.c' files.
// NOTE: we assume the upper and lower inidices defined in those above
// files won't be changing, otherwise we need to change them in Elliptica 
// as well. */
static char *required_field[] = 
{
"alpha",/* lapse: alpha */
"betax","betay","betaz",/* shift: beta^i */
"adm_gxx","adm_gxy","adm_gxz",/* metric: g_ij */
"adm_gyy","adm_gyz","adm_gzz",/* metric: g_ij */
"adm_Kxx","adm_Kxy","adm_Kxz",/* extrinsic curvature: K_ij */
"adm_Kyy","adm_Kyz","adm_Kzz",/* extrinsic curvature: K_ij */
"grhd_vx","grhd_vy","grhd_vz",/* primitive v, from BNSdata v^i defined: 
                              // v^i = u^i/(alpha u^0) + beta^i / alpha */
                         
"grhd_rho",/* primitive rho */
"grhd_p",/* primitive p */
"grhd_epsl",/* primitive epsilon: total_energy_density = grhd_rho(1+grhd_epsl)*/
     
//"adm_rho",/* ADM rho */
//"adm_Sx","adm_Sy","adm_Sz",/* ADM S^i */
//"adm_SSxx","adm_SSxy","adm_SSxz",/* ADM S_ij */
//"adm_SSyy","adm_SSyz","adm_SSzz",/* ADM S_ij */
//"adm_dpsiopsix","adm_dpsiopsiy","adm_dpsiopsiz",
//"adm_psi",
//"adm_ddpsiopsixx","adm_ddpsiopsixy","adm_ddpsiopsixz",
//"adm_ddpsiopsiyy","adm_ddpsiopsiyz","adm_ddpsiopsizz",

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
  
  /* files path */
  sprintf(coords_file_path, "%s/coords_level%d_proc%d.dat", 
    outdir, level->l, rank);
  
  sprintf(fields_file_path, "%s/fields_level%d_proc%d.dat", 
    outdir, level->l, rank);
  
  /* write cartesian coords into file */
  write_coords(level,coords_file_path);
  
  /* call elliptica and interpolate at (x,y,z) given by the coords file */
  call_elliptica_and_write_fields(level,coords_file_path,fields_file_path);
  
  /* populate fields for bam */
  populate_fields_for_bam(level,fields_file_path);
  
  /* delete files */
  if (DELETE)
  {
     char command[STR_LEN_MAX2x]={'\0'};
     int ret;
     
     sprintf(command,"rm -rf %s",coords_file_path);
     printf("System call:\n%s\n",command);
     fflush(stdout);
     ret = system(command);
     printf("System call returned: %d\n",ret);
     fflush(stdout);
     
     sprintf(command,"rm -rf %s",fields_file_path);
     printf("System call:\n%s\n",command);
     fflush(stdout);
     ret = system(command);
     printf("System call returned: %d\n",ret);
     fflush(stdout);
  }
  printf("} Importing initial data from Elliptica --> Done.\n");
  fflush(stdout);
  return 0;
}

/* populating the fields made by Elliptica into bam */
static void populate_fields_for_bam(tL *const level, char *const fields_file_path)
{
  printf("~> Populating BAM variables based on initial data ...\n");
  fflush(stdout);
  FILE *file = 0;
  double *v  = 0;
  char *v_name = 0;
  const double *const x = level->v[Ind("x")];
  const double *const y = level->v[Ind("y")];
  const double *const z = level->v[Ind("z")];
  char *match_str;
  char msg[STR_LEN_MAX];
  int msg_len = (int)strlen(END_MSG)+1;
  int i,f;
  
  /* read fields from the file made by Elliptica */
  file = fopen(fields_file_path,"r");
  if(!file) 
    errorexits("could not open %s", fields_file_path);
    
  /* check if data file is completed */
  fseek(file,-msg_len,SEEK_END);
  ASSERT(fread(msg,(unsigned)msg_len,1,file));
  if (!strstr(msg,END_MSG))
    errorexit("Fields data file is corrupted.\n");
    
  fseek(file,0,SEEK_SET);
  fgets(msg,STR_LEN_MAX,file);/* read first comment */
  
  FReadP_bin(match_str);/* read header */
  if (strcmp(match_str,HEADER))
    errorexit("It could not find the header.\n");
  free(match_str);
  f = 0;
  while(required_field[f])
  {
    /* read field name */
    FReadP_bin(v_name);
    if(strcmp(v_name,required_field[f]))
       errorexit("It could not find the field.\n");
       
    /* free */
    free(v_name);
    
    /* enable */   
    int comp = Ind(required_field[f]);
    enablevarcomp(level, comp);
    v = level->v[comp];
    
    /* read field value */
    forallpoints(level,i)
    {
      FReadV_bin(v[i]);
      ASSERT(isfinite(v[i]));
    }
    
    ++f;
  }
  FReadP_bin(match_str);/* read footer */
  if (strcmp(match_str,FOOTER))
    errorexit("It could not find the footer.\n");
  free(match_str);
  fclose(file);
  
  /* special treatments: */
  /* set psi and its derivs */
  i = Ind("adm_psi"); enablevar(level, i);
  double *psi = level->v[i++];
  forallpoints(level,i)
  {
    psi[i] = 1.0;
  }
  i = Ind("adm_dpsiopsix");   enablevar(level, i);
  i = Ind("adm_ddpsiopsixx"); enablevar(level, i);
  
  /* set S_i's and S_ij's */
  i = Ind("adm_Sx");   enablevar(level, i);
  i = Ind("adm_SSxx"); enablevar(level, i);
  
  /* amd_rho */
  i = Ind("adm_rho"); enablevar(level, i);
  
  double *gxx = Ptr(level, "adm_gxx");
  double *gxy = Ptr(level, "adm_gxy");
  double *gxz = Ptr(level, "adm_gxz");
  double *gyy = Ptr(level, "adm_gyy");
  double *gyz = Ptr(level, "adm_gyz");
  double *gzz = Ptr(level, "adm_gzz");
  
  /* set tr(K0)? */
  if (0)
  {
    i = Ind("adm_K0"); enablevar(level, i);
    double *Kxx = Ptr(level, "adm_Kxx");
    double *Kxy = Ptr(level, "adm_Kxy");
    double *Kxz = Ptr(level, "adm_Kxz");
    double *Kyy = Ptr(level, "adm_Kyy");
    double *Kyz = Ptr(level, "adm_Kyz");
    double *Kzz = Ptr(level, "adm_Kzz");
    double *K0  = Ptr(level, "adm_K0");
    
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
             detg, i, x[i], y[i], z[i]);
    }
  }

  printf("~> Populating BAM variables based on initial data --> Done.\n");
  fflush(stdout);
}

/* call elliptica and interpolate at (x,y,z) given by the file */
static void call_elliptica_and_write_fields(tL *const level,char *const coords_file_path, char *const fields_file_path)
{
  const int rank = bampi_rank();
  const char *const outdir = Gets("outdir");
  FILE *id_parfile = 0;
  char command[STR_LEN_MAX] = {'\0'};
  char elliptica_parfile[STR_LEN_MAX] = {'\0'};
  char id_parfile_path[STR_LEN_MAX] = {'\0'};
  char id_parfile_name[STR_LEN_MAX] = {'\0'};
  char id_outdir[STR_LEN_MAX] = {'\0'};
  char str[STR_LEN_MAX] = {'\0'};
  char *stem;
  int ret,i;
  
  /* initiliazing names and paths */
  sprintf(str,"%s",Gets("EIDdataReader_datadir"));/* --> a/b/c */
  if (str[strlen(str)-1] == '/')/* if there is / at the end of path */
    str[strlen(str)-1] = '\0';/* if a/b/c/ */
  sprintf(id_outdir,"%s",str);
  stem = strrchr(str,'/');/* stem = /c */
  if (!stem)/* if stem = a */
    stem = str;/* stem = a */
  else
    stem++;/* stem = c */
  
  /* copying parfile for this process id */
  sprintf(id_parfile_path,"%s/%s.par",
          id_outdir,stem);/* a/b/c/idfile.par */
  
  /* check if the parfile exists or the name is correct */        
  id_parfile = fopen(id_parfile_path,"r");
  if (!id_parfile)
    errorexits("could not open '%s'.\n "
               "Make sure directory and parfile have the same name.", id_parfile_path);
  fclose(id_parfile);
  
  sprintf(id_parfile_name,"%s_elliptica_level%d_proc%d.par",
          stem,level->l,rank);/* idfile_elliptica_level?_proc?.par */
  sprintf(command,"cp %s %s/%s",
          id_parfile_path,outdir,
          id_parfile_name);/* cp idfile.par idfile_elliptica_level?_proc?.par */
  sprintf(id_parfile_path,"%s/%s",
          outdir,id_parfile_name);/* outdir/idfile_elliptica_level?_proc?.par */
  printf("System call:\n%s\n",command);
  fflush(stdout);
  ret = system(command);
  printf("System call returned: %d\n",ret);
  fflush(stdout);
  
  /* modifying parfile for bam */        
  i = 0;
  str[0] = '\0';
  while(required_field[i])
  {
    strcat(str,required_field[i]);
    strcat(str,",");
    i++;
  }
  id_parfile = fopen(id_parfile_path,"a");
  if (!id_parfile)
    errorexits("could not open %s", id_parfile_path);
    
  /* adding these paramters to Elliptica parfile to direct */
  if(Getv("EIDdateReader_physics",BHNS))
  {
    fprintf(id_parfile,"\n");
    fprintf(id_parfile,BHNS_ EVO_"export_id    = yes\n");
    fprintf(id_parfile,Pbhns_"export_id        = yes\n");
    fprintf(id_parfile,Pbhns_"coords_file_path = %s\n",coords_file_path);
    fprintf(id_parfile,Pbhns_"fields_file_path = %s\n",fields_file_path);
    fprintf(id_parfile,Pbhns_"fields_name      = %s\n",str);
    fprintf(id_parfile,Pbhns_"filler_method    = %s\n",Gets("EIDdateReader_BHfiller"));
    fprintf(id_parfile,Pbhns_"checkpoint_file_path     = %s/checkpoint.dat\n",id_outdir);
    fprintf(id_parfile,BHNS_ EVO_"checkpoint_file_path = %s/checkpoint.dat\n",id_outdir);
    fprintf(id_parfile,"\n");
  }
  else if(Getv("EIDdateReader_physics",SBH))
  {
    fprintf(id_parfile,"\n");
    fprintf(id_parfile,SBH_ EVO_"export_id    = yes\n");
    fprintf(id_parfile,Psbh_"export_id        = yes\n");
    fprintf(id_parfile,Psbh_"coords_file_path = %s\n",coords_file_path);
    fprintf(id_parfile,Psbh_"fields_file_path = %s\n",fields_file_path);
    fprintf(id_parfile,Psbh_"fields_name      = %s\n",str);
    fprintf(id_parfile,Psbh_"filler_method    = %s\n",Gets("EIDdateReader_BHfiller"));
    fprintf(id_parfile,Psbh_"checkpoint_file_path     = %s/checkpoint.dat\n",id_outdir);
    fprintf(id_parfile,SBH_ EVO_"checkpoint_file_path = %s/checkpoint.dat\n",id_outdir);
    fprintf(id_parfile,"\n");
  }
  else
    errorexits("No such %s implemented!",Gets("EIDdateReader_physics"));

  fclose(id_parfile);
  
  /* call Elliptica */
  sprintf(command, "%s %s",Gets("EIDdataReader_exe"),id_parfile_path);
  printf("System call:\n%s\n", command);
  fflush(stdout);
  ret = system(command);
  printf("System call returned: %d\n", ret);
  fflush(stdout);

  /* remove par id_parfile */
  if (DELETE)
  {
    /* rm -rf outdir/idfile_elliptica_level?_proc?.par */
    sprintf(command,"rm -rf %s",id_parfile_path);
    printf("System call:\n%s\n",command);
    fflush(stdout);
    ret = system(command);
  }
  printf("System call returned: %d\n", ret);
  fflush(stdout);
}
  
/* write grid point in Cartesian coords. (x,y,z) into a file.
// ->return value: pointer to the file path */
static void write_coords(tL *const level,char *const coords_file_path)
{
  FILE *file = 0;
  char title_line[STR_LEN_MAX];
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  const double *const x = level->v[Ind("x")];
  const double *const y = level->v[Ind("y")];
  const double *const z = level->v[Ind("z")];
  const double x_CM = Getd("EIDdateReader_x_CM");
  const double y_CM = Getd("EIDdateReader_y_CM");
  const double z_CM = Getd("EIDdateReader_z_CM");
  int i;
  
  file = fopen(coords_file_path, "wb");
  if(!file) 
    errorexits("could not open %s", coords_file_path);
  
  /* write coords file */
  fprintf(file,"# this file contains (x,y,z) of bam grid points.\n");
  sprintf(title_line,"%s",HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  FWriteV_bin(level->npoints,1);
  /* note: r_elliptica = r_CM + r_bam = > r_bam = r_elliptica - r_CM. */
  forallpoints(level, i)
  {
    double x_ell = x[i]+x_CM;
    double y_ell = y[i]+y_CM;
    double z_ell = z[i]+z_CM;
    
    FWriteV_bin(x_ell,1);
    FWriteV_bin(y_ell,1);
    FWriteV_bin(z_ell,1);
  }
  sprintf(title_line,"%s",FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  fclose(file);
}

/* set some pars important for making the grid */
int EIDpreGrid(tL *const level)
{
  /* set pre grid parameters */
  printf("Setting some parameters relevant for setting up bam's grid:\n");
  
  if(Getv("EIDdateReader_physics",BHNS))
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
    Setd("EIDdateReader_x_CM",x_CM);
    Setd("EIDdateReader_y_CM",y_CM);
    Setd("EIDdateReader_z_CM",z_CM);
    
    READ_PARAMETER_FROM_FILE(ns_center_x,"NS_center_x")
    READ_PARAMETER_FROM_FILE(ns_center_y,"NS_center_y")
    READ_PARAMETER_FROM_FILE(ns_center_z,"NS_center_z")
    
    READ_PARAMETER_FROM_FILE(bh_center_x,"BH_center_x")
    READ_PARAMETER_FROM_FILE(bh_center_y,"BH_center_y")
    READ_PARAMETER_FROM_FILE(bh_center_z,"BH_center_z")
    
    READ_PARAMETER_FROM_FILE(ns_m,"NS_baryonic_mass")
    READ_PARAMETER_FROM_FILE(bh_m,"BH_irreducible_mass")
    
    /* set pre grid params 
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    Setd("mass1", ns_m);
    Setd("mass2", bh_m);
    Setd("px1", ns_center_x-x_CM);
    Setd("py1", ns_center_y-y_CM);
    Setd("pz1", ns_center_z-z_CM);
    Setd("px2", bh_center_x-x_CM);
    Setd("py2", bh_center_y-y_CM);
    Setd("pz2", bh_center_z-z_CM);
    
    fclose(file);
  }
  else if(Getv("EIDdateReader_physics",SBH))
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
    Setd("EIDdateReader_x_CM",x_CM);
    Setd("EIDdateReader_y_CM",y_CM);
    Setd("EIDdateReader_z_CM",z_CM);
    
    READ_PARAMETER_FROM_FILE(bh_center_x,"BH_center_x")
    READ_PARAMETER_FROM_FILE(bh_center_y,"BH_center_y")
    READ_PARAMETER_FROM_FILE(bh_center_z,"BH_center_z")
    
    READ_PARAMETER_FROM_FILE(bh_m,"BH_irreducible_mass")
    
    /* set pre grid params 
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    Setd("mass1", bh_m);
    Setd("mass2", 0);
    Setd("px1", bh_center_x-x_CM);
    Setd("py1", bh_center_y-y_CM);
    Setd("pz1", bh_center_z-z_CM);
    
    fclose(file);
  }
  else
    errorexits("No such %s implemented!",Gets("EIDdateReader_physics"));
  
  return 0;
}
