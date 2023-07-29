/* Alireza Rashti Feb 2020 */

#include "bam.h"
#include "elliptica_id_reader_lib.h"
#include "EIDdataReader.h"

static void read_fields_from_file(tL *const level, char *const fields_file_path);
static void interpolate_field(tL *const level);
static void debug_save_fields_to_txt(tL *const level,const char *const field_names[]);

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
// 1. write the grid points in Cartesian coords., (x,y,z), into 1D arrays
// 2. specify the required fields to be interpolated for bam
// 3. call Elliptica to interpolate the fields on these (x,y,z)'s
// 4. read the interpolated values and populate Bam's fields
//
// NOTE: the format of file names containing fields are "fields_level%d_proc%d.dat"
//
// ->return value: 0 */
int EIDdataReader(tL *const level)
{
  printf("{ Importing initial data from Elliptica ...\n");

  const char *const outdir = Gets("outdir");
  const int rank = bampi_rank();
  char fields_file_path[STR_LEN_MAX] = {'\0'};

  // if this is a checkpoint, just read field values from the save files
  if(Getv("checkpoint","restart") && Getv("EIDdataReader_save","yes"))
  {
    sprintf(fields_file_path, "%s_previous/%s/fields_level%d_proc%d.dat",
      outdir,Gets("EIDdataReader_outdir"), level->l, rank);

    /* populate fields for bam */
    read_fields_from_file(level,fields_file_path);
  }
  // load from the specified directory
  else if (!Getv("EIDdataReader_loadfrom","not_set!"))
  {
    sprintf(fields_file_path, "%s/fields_level%d_proc%d.dat",
        Gets("EIDdataReader_loadfrom"), level->l, rank);
    
    /* populate fields for bam */
    read_fields_from_file(level,fields_file_path);
  }
  // interpolate from scratch and save (if asked)
  else
  {
    interpolate_field(level);
  }
  
  printf("} Importing initial data from Elliptica --> Done.\n");
  fflush(stdout);
  return 0;
}

/* set some params that important for making the grid */
int EIDpreGrid(tL *const level)
{
  /* set pre grid parameters */
  printf("Setting some parameters relevant for setting up bam's grid:\n");
  
  /* initialize ID Reader */ 
  Elliptica_ID_Reader_T *idr = 
    elliptica_id_reader_init(Gets("EIDdataReader_checkpoint"),"generic");
  
  if(Getv("EIDdataReader_physics","BHNS"))
  {
    const double sys_cm[3] = {
      idr->get_param_dbl("BHNS_x_CM",idr),
      idr->get_param_dbl("BHNS_y_CM",idr),
      idr->get_param_dbl("BHNS_z_CM",idr)
    };
    
    const double ns_c[3] = {
      idr->get_param_dbl("NS_center_x",idr),
      idr->get_param_dbl("NS_center_y",idr),
      idr->get_param_dbl("NS_center_z",idr)
    };

    const double bh_c[3] = {
      idr->get_param_dbl("BH_center_x",idr),
      idr->get_param_dbl("BH_center_y",idr),
      idr->get_param_dbl("BH_center_z",idr)
    };

    const double ns_m = idr->get_param_dbl("NS_baryonic_mass_current",idr);
    const double bh_m = idr->get_param_dbl("BH_irreducible_mass_current",idr);
    
    // set pre grid params:

    MySetd("mass1", ns_m);
    MySetd("mass2", bh_m);
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    MySetd("px1", ns_c[0]-sys_cm[0]);
    MySetd("py1", ns_c[1]-sys_cm[1]);
    MySetd("pz1", ns_c[2]-sys_cm[2]);
    MySetd("px2", bh_c[0]-sys_cm[0]);
    MySetd("py2", bh_c[1]-sys_cm[1]);
    MySetd("pz2", bh_c[2]-sys_cm[2]);
    
    /* for the bam BHfiller */
    if (1)
    {
      /* note: the notation different here and object 1 is bh. */
      MySetd("bhmass1",bh_m);
      MySetd("bhmass2",0.);
      MySetd("bhx1", Getd("px2"));
      MySetd("bhy1", Getd("py2"));
      MySetd("bhz1", Getd("pz2"));
    }
  }
  else if(Getv("EIDdataReader_physics","SBH"))
  {
    const double sys_cm[3] = {
      idr->get_param_dbl("SBH_x_CM",idr),
      idr->get_param_dbl("SBH_y_CM",idr),
      idr->get_param_dbl("SBH_z_CM",idr)
    };
    
    const double bh_c[3] = {
      idr->get_param_dbl("BH_center_x",idr),
      idr->get_param_dbl("BH_center_y",idr),
      idr->get_param_dbl("BH_center_z",idr)
    };

    const double bh_m = idr->get_param_dbl("BH_irreducible_mass_current",idr);
    
    // set pre grid params:

    MySetd("mass1", bh_m);
    MySetd("mass2", 0.);
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    MySetd("px1", bh_c[0]-sys_cm[0]);
    MySetd("py1", bh_c[1]-sys_cm[1]);
    MySetd("pz1", bh_c[2]-sys_cm[2]);
    
    /* for the bam BHfiller */
    if (1)
    {
      MySetd("bhmass1",bh_m);
      MySetd("bhmass2",0.);
      MySetd("bhx1", Getd("px1"));
      MySetd("bhy1", Getd("py1"));
      MySetd("bhz1", Getd("pz1"));
    }
  }
  else if(Getv("EIDdataReader_physics","NSNS"))
  {
    const double sys_cm[3] = {
      idr->get_param_dbl("NSNS_x_CM",idr),
      idr->get_param_dbl("NSNS_y_CM",idr),
      idr->get_param_dbl("NSNS_z_CM",idr)
    };
    
    const double ns1_c[3] = {
      idr->get_param_dbl("NS1_center_x",idr),
      idr->get_param_dbl("NS1_center_y",idr),
      idr->get_param_dbl("NS1_center_z",idr)
    };

    const double ns2_c[3] = {
      idr->get_param_dbl("NS2_center_x",idr),
      idr->get_param_dbl("NS2_center_y",idr),
      idr->get_param_dbl("NS2_center_z",idr)
    };

    const double ns1_m = idr->get_param_dbl("NS1_baryonic_mass_current",idr);
    const double ns2_m = idr->get_param_dbl("NS2_baryonic_mass_current",idr);
    
    // set pre grid params:

    MySetd("mass1", ns1_m);
    MySetd("mass2", ns2_m);
    // note: r_elliptica = r_CM + r_bam 
    //       = > r_bam = r_elliptica - r_CM. */
    MySetd("px1", ns1_c[0]-sys_cm[0]);
    MySetd("py1", ns1_c[1]-sys_cm[1]);
    MySetd("pz1", ns1_c[2]-sys_cm[2]);
    MySetd("px2", ns2_c[0]-sys_cm[0]);
    MySetd("py2", ns2_c[1]-sys_cm[1]);
    MySetd("pz2", ns2_c[2]-sys_cm[2]);
    
  }
  else
  {
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  }
   
  // free workspace
  elliptica_id_reader_free(idr);
  
  return 0;
}

/* populating the fields made by Elliptica into bam */
static void read_fields_from_file(tL *const level, char *const fields_file_path)
{
  printf("calling: %s ...\n",__func__);
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
  
  /* read fields content */
  if(Getv("EIDdataReader_physics","BHNS") || Getv("EIDdataReader_physics","NSNS"))
  {
    f = 0;
    while(import_fields_with_matter[f])
    {
      /* read field name */
      FReadP_bin(v_name);
      if(strcmp(v_name,import_fields_with_matter[f]))
         errorexit("It could not find the field.\n");
         
      /* free */
      free(v_name);
      
      /* enable */   
      int comp = Ind(import_fields_with_matter[f]);
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
  }
  else if(Getv("EIDdataReader_physics","SBH"))
  {
    f = 0;
    while(import_fields_no_matter[f])
    {
      /* read field name */
      FReadP_bin(v_name);
      if(strcmp(v_name,import_fields_no_matter[f]))
         errorexit("It could not find the field.\n");
         
      /* free */
      free(v_name);
      
      /* enable */   
      int comp = Ind(import_fields_no_matter[f]);
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
  }
  else
  {
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
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
  
  if(Getv("EIDdataReader_physics","BHNS"))
  {
  /* set S_i's and S_ij's */
    i = Ind("adm_Sx");   enablevar(level, i);
    i = Ind("adm_SSxx"); enablevar(level, i);
    
    /* adm_rho */
    i = Ind("adm_rho"); enablevar(level, i);
  }
  
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

}

// after interpolation of *field_names* save the data into the given *file_path*
static void save_interpolated_values(tL *const level,const char *const file_path, 
                                     const char *const field_names[])
{
  printf("calling: %s ...\n",__func__);
  fflush(stdout);

  FILE *file = 0;
  char title_line[STR_LEN_MAX];
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  char msg[STR_LEN_MAX] = {'\0'};
  char *const p_msg = msg;/* to avoid GCC warning for FWriteP_bin */
  char fnames[STR_LEN_MAX] = {'\0'};
  int fld,i;
  
  // open
  file = fopen(file_path,"wb");
  if(!file) errorexits("could not open %s", file_path);
  
  // prepare and write the header
  fld = 0;
  fnames[0] = '\0';
  while(field_names[fld])
  {
    strcat(fnames,field_names[fld]);
    strcat(fnames,",");
    fld++;
  }
  fnames[strlen(fnames)-1] = '\0';/* rm the last ',' */
  fprintf(file,"# this file contains values of %s\n",fnames);
  sprintf(title_line,"%s",HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);

  // write
  fld = 0;
  while(field_names[fld])
  {
    // write field name
    FWriteP_bin(field_names[fld],strlen(field_names[fld])+1);
    const double *interp_v = level->v[Ind(field_names[fld])];
    
    // write its values
    forallpoints(level, i)
    {
      /* write it into the fields_file */
      FWriteV_bin(interp_v[i],1);
    }
    fld++;
  }
  
  // close 
  sprintf(title_line,"%s",FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  sprintf(msg,"%s",END_MSG);
  FWriteP_bin(p_msg,strlen(msg)+1);
  Fclose(file);
}

// interpolate from scratch and save (if asked)
static void interpolate_field(tL *const level)
{
  printf("calling: %s ...\n",__func__);
  fflush(stdout);

  double *const x = level->v[Ind("x")];
  double *const y = level->v[Ind("y")];
  double *const z = level->v[Ind("z")];
  const int rank  = bampi_rank();
  const char *const outdir = Gets("outdir");
  char fields_file_path[STR_LEN_MAX] = {'\0'};
  int i;
  
  /* initialize ID Reader */ 
  Elliptica_ID_Reader_T *idr = 
    elliptica_id_reader_init(Gets("EIDdataReader_checkpoint"),"generic");
  
  // specify interpolation points
  idr->npoints  = level->npoints;
  idr->x_coords = x;
  idr->y_coords = y;
  idr->z_coords = z;
  
  /* specify interpolating fields */
  if(Getv("EIDdataReader_physics","BHNS"))
  {
    idr->ifields = "alpha,betax,betay,betaz,"
                   "adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,"
                   "adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,"
                   "grhd_rho,grhd_p,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
    idr->set_param("BHNS_filler_method",Gets("EIDdataReader_BHfiller"),idr);
    idr->set_param("ADM_B1I_form","zero",idr); 
  }
  else if(Getv("EIDdataReader_physics","NSNS"))
  {
    idr->ifields = "alpha,betax,betay,betaz,"
                   "adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,"
                   "adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,"
                   "grhd_rho,grhd_p,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
    idr->set_param("ADM_B1I_form","zero",idr); 
  }
  else if(Getv("EIDdataReader_physics","SBH"))
  {
    idr->ifields = "alpha,betax,betay,betaz,"
                   "adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,"
                   "adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz";
  }
  else
  {
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  }
  
  /* interpolate */
  elliptica_id_reader_interpolate(idr);
  
  // save level indices for a slight optimization (no matter field)
  // bam's indices
  const int ibam_alpha = Ind("alpha");
  const int ibam_betax = Ind("betax");
  const int ibam_betay = Ind("betay");
  const int ibam_betaz = Ind("betaz");

  const int ibam_adm_gxx = Ind("adm_gxx");
  const int ibam_adm_gxy = Ind("adm_gxy");
  const int ibam_adm_gxz = Ind("adm_gxz");
  const int ibam_adm_gyy = Ind("adm_gyy");
  const int ibam_adm_gyz = Ind("adm_gyz");
  const int ibam_adm_gzz = Ind("adm_gzz");
  
  const int ibam_adm_Kxx = Ind("adm_Kxx");
  const int ibam_adm_Kxy = Ind("adm_Kxy");
  const int ibam_adm_Kxz = Ind("adm_Kxz");
  const int ibam_adm_Kyy = Ind("adm_Kyy");
  const int ibam_adm_Kyz = Ind("adm_Kyz");
  const int ibam_adm_Kzz = Ind("adm_Kzz");

  // elliptica's indices
  const int iell_alpha   = idr->indx("alpha");
  const int iell_betax   = idr->indx("betax");
  const int iell_betay   = idr->indx("betay");
  const int iell_betaz   = idr->indx("betaz");
  
  const int iell_adm_gxx = idr->indx("adm_gxx");
  const int iell_adm_gxy = idr->indx("adm_gxy");
  const int iell_adm_gxz = idr->indx("adm_gxz");
  const int iell_adm_gyy = idr->indx("adm_gyy");
  const int iell_adm_gyz = idr->indx("adm_gyz");
  const int iell_adm_gzz = idr->indx("adm_gzz");

  const int iell_adm_Kxx = idr->indx("adm_Kxx");
  const int iell_adm_Kxy = idr->indx("adm_Kxy");
  const int iell_adm_Kxz = idr->indx("adm_Kxz");
  const int iell_adm_Kyy = idr->indx("adm_Kyy");
  const int iell_adm_Kyz = idr->indx("adm_Kyz");
  const int iell_adm_Kzz = idr->indx("adm_Kzz");
  
  // populate fields
  if(Getv("EIDdataReader_physics","BHNS") || Getv("EIDdataReader_physics","NSNS") )
  {
    // bam's indices
    const int ibam_grhd_rho  = Ind("grhd_rho");
    const int ibam_grhd_epsl = Ind("grhd_epsl");
    const int ibam_grhd_p    = Ind("grhd_p");
    const int ibam_grhd_vx   = Ind("grhd_vx");
    const int ibam_grhd_vy   = Ind("grhd_vy");
    const int ibam_grhd_vz   = Ind("grhd_vz");

    // elliptica's indices
    const int iell_grhd_rho  = idr->indx("grhd_rho");
    const int iell_grhd_epsl = idr->indx("grhd_epsl");
    const int iell_grhd_p    = idr->indx("grhd_p");
    const int iell_grhd_vx   = idr->indx("grhd_vx");
    const int iell_grhd_vy   = idr->indx("grhd_vy");
    const int iell_grhd_vz   = idr->indx("grhd_vz");

    int fld = 0;
    while(import_fields_with_matter[fld])
    {
     /* enable */   
      int comp = Ind(import_fields_with_matter[fld]);
      enablevarcomp(level, comp);
      fld++;
    }

    forallpoints(level, i)
    {
      level->v[ibam_alpha][i]   = idr->field[iell_alpha][i];
      level->v[ibam_betax][i]   = idr->field[iell_betax][i];
      level->v[ibam_betay][i]   = idr->field[iell_betay][i];
      level->v[ibam_betaz][i]   = idr->field[iell_betaz][i];

      level->v[ibam_adm_gxx][i] = idr->field[iell_adm_gxx][i];
      level->v[ibam_adm_gxy][i] = idr->field[iell_adm_gxy][i];
      level->v[ibam_adm_gxz][i] = idr->field[iell_adm_gxz][i];
      level->v[ibam_adm_gyy][i] = idr->field[iell_adm_gyy][i];
      level->v[ibam_adm_gyz][i] = idr->field[iell_adm_gyz][i];
      level->v[ibam_adm_gzz][i] = idr->field[iell_adm_gzz][i];

      level->v[ibam_adm_Kxx][i] = idr->field[iell_adm_Kxx][i];
      level->v[ibam_adm_Kxy][i] = idr->field[iell_adm_Kxy][i];
      level->v[ibam_adm_Kxz][i] = idr->field[iell_adm_Kxz][i];
      level->v[ibam_adm_Kyy][i] = idr->field[iell_adm_Kyy][i];
      level->v[ibam_adm_Kyz][i] = idr->field[iell_adm_Kyz][i];
      level->v[ibam_adm_Kzz][i] = idr->field[iell_adm_Kzz][i];

      level->v[ibam_grhd_rho][i]  = idr->field[iell_grhd_rho][i];
      level->v[ibam_grhd_epsl][i] = idr->field[iell_grhd_epsl][i];
      level->v[ibam_grhd_p][i]    = idr->field[iell_grhd_p][i];
      level->v[ibam_grhd_vx][i]   = idr->field[iell_grhd_vx][i];
      level->v[ibam_grhd_vy][i]   = idr->field[iell_grhd_vy][i];
      level->v[ibam_grhd_vz][i]   = idr->field[iell_grhd_vz][i];
    }
    
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
    /* adm_rho */
    i = Ind("adm_rho"); enablevar(level, i);
    
    if (Getv("EIDdataReader_save","yes"))
    {
      sprintf(fields_file_path, "%s/%s/fields_level%d_proc%d.dat", 
        outdir,Gets("EIDdataReader_outdir"),level->l, rank);
          
      save_interpolated_values(level,fields_file_path,import_fields_with_matter);
    }
    
    if (DEBUG_TXT)
    {
      debug_save_fields_to_txt(level,import_fields_with_matter);
    }
    
  }
  else if(Getv("EIDdataReader_physics","SBH"))
  {
    int fld = 0;
    while(import_fields_no_matter[fld])
    {
     /* enable */   
      int comp = Ind(import_fields_no_matter[fld]);
      enablevarcomp(level, comp);
      fld++;
    }

    forallpoints(level, i)
    {
      level->v[ibam_alpha][i]   = idr->field[iell_alpha][i];
      level->v[ibam_betax][i]   = idr->field[iell_betax][i];
      level->v[ibam_betay][i]   = idr->field[iell_betay][i];
      level->v[ibam_betaz][i]   = idr->field[iell_betaz][i];

      level->v[ibam_adm_gxx][i] = idr->field[iell_adm_gxx][i];
      level->v[ibam_adm_gxy][i] = idr->field[iell_adm_gxy][i];
      level->v[ibam_adm_gxz][i] = idr->field[iell_adm_gxz][i];
      level->v[ibam_adm_gyy][i] = idr->field[iell_adm_gyy][i];
      level->v[ibam_adm_gyz][i] = idr->field[iell_adm_gyz][i];
      level->v[ibam_adm_gzz][i] = idr->field[iell_adm_gzz][i];

      level->v[ibam_adm_Kxx][i] = idr->field[iell_adm_Kxx][i];
      level->v[ibam_adm_Kxy][i] = idr->field[iell_adm_Kxy][i];
      level->v[ibam_adm_Kxz][i] = idr->field[iell_adm_Kxz][i];
      level->v[ibam_adm_Kyy][i] = idr->field[iell_adm_Kyy][i];
      level->v[ibam_adm_Kyz][i] = idr->field[iell_adm_Kyz][i];
      level->v[ibam_adm_Kzz][i] = idr->field[iell_adm_Kzz][i];
    }
    
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

    if (Getv("EIDdataReader_save","yes"))
    {
      sprintf(fields_file_path, "%s/%s/fields_level%d_proc%d.dat", 
        outdir,Gets("EIDdataReader_outdir"),level->l, rank);
      
      save_interpolated_values(level,fields_file_path,import_fields_no_matter);
    }
    
    if (DEBUG_TXT)
    {
      debug_save_fields_to_txt(level,import_fields_no_matter);
    }
  }
  else
  {
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  }

  /* some tests */
  double *const gxx = level->v[ibam_adm_gxx];
  double *const gxy = level->v[ibam_adm_gxy];
  double *const gxz = level->v[ibam_adm_gxz];
  double *const gyy = level->v[ibam_adm_gyy];
  double *const gyz = level->v[ibam_adm_gyz];
  double *const gzz = level->v[ibam_adm_gzz];
  double detg;

  /* check if determinant is positvie, since the metric is Reimannian */
  forallpoints(level,i)
  {
    detg=(2.*gxy[i]*gxz[i]*gyz[i] + gxx[i]*gyy[i]*gzz[i] -
             gzz[i]*gxy[i]*gxy[i] - gyy[i]*gxz[i]*gxz[i] -
             gxx[i]*gyz[i]*gyz[i]);

    if(detg<=0.)
    {
      printf("det(g_ij)=%g at ccc=i=%d:  x=%g y=%g z=%g\n", 
             detg, i, x[i], y[i], z[i]);
    }
  }
  // free workspace
  elliptica_id_reader_free(idr);
}

// save field to txt files for debug purposes
static void debug_save_fields_to_txt(tL *const level,const char *const field_names[])
{
  printf("calling: %s ...\n",__func__);
  fflush(stdout);
  
  const char *const outdir = Gets("outdir");
  const int rank = bampi_rank();
  FILE *file = 0;
  char file_path[STR_LEN_MAX] = {'\0'};
  int fld,i;
  
  fld = 0;
  while(field_names[fld])
  {
    sprintf(file_path, "%s/%s/debug_%s_level%d_proc%d.txt", 
            outdir,Gets("EIDdataReader_outdir"),field_names[fld],level->l, rank);
  
    file = fopen(file_path,"w");
    if(!file) errorexits("could not open %s", file_path);
  
    const double *interp_v = level->v[Ind(field_names[fld])];
    fprintf(file,"# point_index %s\n",field_names[fld]);
    // write its values
    forallpoints(level, i)
    {
      fprintf(file,"%d %0.15e\n",i,interp_v[i]);
    }
    fld++;
    
    Fclose(file);
  }
}
