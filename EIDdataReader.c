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
// 1. write the grid points in Cartesian coords., (x,y,z), into 1D arrays
// 2. specify the required fields to be interpolated for bam
// 3. call Elliptica to interpolate the fields on these (x,y,z)'s
// 4. read the interpolated values and populate Bam's fields
// ->return value: 0 */
int EIDdataReader(tL *const level)
{
  printf("{ Importing initial data from Elliptica ...\n");
  
  const char *checkpnt  = Gets("EIDdataReader_checkpoint");

  /* initialize ID Reader */ 
  Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(checkpnt_path,"generic");
  idr->npoints  = level->npoints;
  idr->x_coords = level->v[Ind("x")];
  idr->y_coords = level->v[Ind("y")];
  idr->z_coords = level->v[Ind("z")];

  /* read fields content */
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
  
  printf("~> Populating BAM variables based on initial data ...\n");
  fflush(stdout);
  
  if(Getv("EIDdataReader_physics","BHNS") || Getv("EIDdataReader_physics","NSNS") )
  {
    const int ialpha = idr->indx("alpha");
    const int ibetax = idr->indx("betax");
    const int ibetay = idr->indx("betay");
    const int ibetaz = idr->indx("betaz");

    const int iadm_gxx   = idr->indx("adm_gxx");
    const int iadm_gxy   = idr->indx("adm_gxy");
    const int iadm_gxz   = idr->indx("adm_gxz");
    const int iadm_gyy   = idr->indx("adm_gyy");
    const int iadm_gyz   = idr->indx("adm_gyz");
    const int iadm_gzz   = idr->indx("adm_gzz");

    const int iadm_Kxx   = idr->indx("adm_Kxx");
    const int iadm_Kxy   = idr->indx("adm_Kxy");
    const int iadm_Kxz   = idr->indx("adm_Kxz");
    const int iadm_Kyy   = idr->indx("adm_Kyy");
    const int iadm_Kyz   = idr->indx("adm_Kyz");
    const int iadm_Kzz   = idr->indx("adm_Kzz");

    const int igrhd_rho   = idr->indx("grhd_rho");
    const int igrhd_epsl  = idr->indx("grhd_epsl");
    const int igrhd_p     = idr->indx("grhd_p");
    const int igrhd_vx    = idr->indx("grhd_vx");
    const int igrhd_vy    = idr->indx("grhd_vy");
    const int igrhd_vz    = idr->indx("grhd_vz");

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
      Ptr(level, "alpha")[i]   = idr->field[ialpha][i];
      Ptr(level, "betax")[i]   = idr->field[ibetax][i];
      Ptr(level, "betay")[i]   = idr->field[ibetay][i];
      Ptr(level, "betaz")[i]   = idr->field[ibetaz][i];

      Ptr(level, "adm_gxx")[i] = idr->field[iadm_gxx][i];
      Ptr(level, "adm_gxy")[i] = idr->field[iadm_gxy][i];
      Ptr(level, "adm_gxz")[i] = idr->field[iadm_gxz][i];
      Ptr(level, "adm_gyy")[i] = idr->field[iadm_gyy][i];
      Ptr(level, "adm_gyz")[i] = idr->field[iadm_gyz][i];
      Ptr(level, "adm_gzz")[i] = idr->field[iadm_gzz][i];

      Ptr(level, "adm_Kxx")[i] = idr->field[iadm_Kxx][i];
      Ptr(level, "adm_Kxy")[i] = idr->field[iadm_Kxy][i];
      Ptr(level, "adm_Kxz")[i] = idr->field[iadm_Kxz][i];
      Ptr(level, "adm_Kyy")[i] = idr->field[iadm_Kyy][i];
      Ptr(level, "adm_Kyz")[i] = idr->field[iadm_Kyz][i];
      Ptr(level, "adm_Kzz")[i] = idr->field[iadm_Kzz][i];

      Ptr(level, "grhd_rho")[i]  = idr->field[igrhd_rho][i];
      Ptr(level, "grhd_epsl")[i] = idr->field[igrhd_epsl][i];
      Ptr(level, "grhd_p")[i]    = idr->field[igrhd_p][i];
      Ptr(level, "grhd_vx")[i]   = idr->field[igrhd_vx][i];
      Ptr(level, "grhd_vy")[i]   = idr->field[igrhd_vy][i];
      Ptr(level, "grhd_vz")[i]   = idr->field[igrhd_vz][i];
    }
  }
  else if(Getv("EIDdataReader_physics","SBH"))
  {
    const int ialpha = idr->indx("alpha");
    const int ibetax = idr->indx("betax");
    const int ibetay = idr->indx("betay");
    const int ibetaz = idr->indx("betaz");

    const int iadm_gxx   = idr->indx("adm_gxx");
    const int iadm_gxy   = idr->indx("adm_gxy");
    const int iadm_gxz   = idr->indx("adm_gxz");
    const int iadm_gyy   = idr->indx("adm_gyy");
    const int iadm_gyz   = idr->indx("adm_gyz");
    const int iadm_gzz   = idr->indx("adm_gzz");

    const int iadm_Kxx   = idr->indx("adm_Kxx");
    const int iadm_Kxy   = idr->indx("adm_Kxy");
    const int iadm_Kxz   = idr->indx("adm_Kxz");
    const int iadm_Kyy   = idr->indx("adm_Kyy");
    const int iadm_Kyz   = idr->indx("adm_Kyz");
    const int iadm_Kzz   = idr->indx("adm_Kzz");

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
      Ptr(level, "alpha")[i]   = idr->field[ialpha][i];
      Ptr(level, "betax")[i]   = idr->field[ibetax][i];
      Ptr(level, "betay")[i]   = idr->field[ibetay][i];
      Ptr(level, "betaz")[i]   = idr->field[ibetaz][i];

      Ptr(level, "adm_gxx")[i] = idr->field[iadm_gxx][i];
      Ptr(level, "adm_gxy")[i] = idr->field[iadm_gxy][i];
      Ptr(level, "adm_gxz")[i] = idr->field[iadm_gxz][i];
      Ptr(level, "adm_gyy")[i] = idr->field[iadm_gyy][i];
      Ptr(level, "adm_gyz")[i] = idr->field[iadm_gyz][i];
      Ptr(level, "adm_gzz")[i] = idr->field[iadm_gzz][i];

      Ptr(level, "adm_Kxx")[i] = idr->field[iadm_Kxx][i];
      Ptr(level, "adm_Kxy")[i] = idr->field[iadm_Kxy][i];
      Ptr(level, "adm_Kxz")[i] = idr->field[iadm_Kxz][i];
      Ptr(level, "adm_Kyy")[i] = idr->field[iadm_Kyy][i];
      Ptr(level, "adm_Kyz")[i] = idr->field[iadm_Kyz][i];
      Ptr(level, "adm_Kzz")[i] = idr->field[iadm_Kzz][i];
    }
  }
  else
  {
    errorexits("No such %s implemented!",Gets("EIDdataReader_physics"));
  }

  printf("~> Populating BAM variables based on initial data --> Done.\n");
  fflush(stdout);
  
  // free workspace
  elliptica_id_reader_free(idr);
  printf("} Importing initial data from Elliptica --> Done.\n");
  fflush(stdout);

  return 0;
}
