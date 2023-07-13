#include <assert.h>

/* bothered with compiler warning for const qualifier? */
#define const 

/* file and parameter prefix, let's keep them capitalized for Elliptica */
#define BHNS_ "BHNS_"
#define SBH_  "SBH_"
#define BHNS  "BHNS"
#define SBH   "SBH"
#define EVO_  "BAM_"

/* instruct Elliptica to save new changes in checkpoint parameter */
#define MODIFY "modify:"

/* handy macro */
#define Pbhns_ MODIFY BHNS_ EVO_
#define Psbh_  MODIFY SBH_  EVO_

#define STR_LEN_MAX   (1000)
#define STR_LEN_MAX2x (2000)

#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"

#define ASSERT assert

/* deleting files created by Elliptica */
#define DELETE (1)

/* this is how we write binary data: first write size and then value. 
// thus, when we wanna read the data the first one gives of the memory allocation 
// and the next gives us value: */

/* write pointer: */
#define FWriteP_bin(x,y) \
if (x){\
  unsigned SIZE_ = (unsigned)(y);\
  ASSERT(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
  ASSERT(fwrite(x,sizeof(*(x)),SIZE_,file));\
}else{\
  unsigned SIZE_ = 0;\
  ASSERT(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
}

/* write variable: */
#define FWriteV_bin(x,y) \
{\
  unsigned SIZE_ = (unsigned)(y);\
  ASSERT(fwrite(&SIZE_,sizeof(SIZE_),1,file));\
  ASSERT(fwrite(&(x),sizeof(x),SIZE_,file));\
}

/* read pointer */
#define FReadP_bin(x) {\
  unsigned SIZE_ = 0;\
  ASSERT(fread(&SIZE_, sizeof(SIZE_),1,file));\
  if (SIZE_) {\
    x = calloc(SIZE_,sizeof(*(x))),ASSERT(x);\
    ASSERT(fread(x,sizeof(*(x)),SIZE_,file));}\
  else { x = 0;}}
  
/* read variable */
#define FReadV_bin(x) {\
  unsigned SIZE_ = 0;\
  ASSERT(fread(&SIZE_, sizeof(SIZE_),1,file));\
  ASSERT(fread(&(x),sizeof(x),SIZE_,file));}

/* get parameter */
#define READ_PARAMETER_FROM_FILE(y,x) \
  fseek(file,0,SEEK_SET);\
  ret = fgetparameter(file, x, str);\
  if (ret == EOF)\
    errorexits("could not find %s parameter", x);\
  y = atof(str);\
  printf("%-30s = %+g\n",x,y);

/* Setd and print */
#define MySetd(x,y) \
  Setd(x,(y)); \
  printf("%-30s = %+g\n",x,(y));

int EIDdataReader(tL *const level);
static void populate_fields_for_bam(tL *const level, char *const checkpoint_dir);





