#include <assert.h>

/* Setd and print */
#define MySetd(x,y) \
  Setd(x,(y)); \
  printf("%-30s = %+g\n",x,(y));

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

#define STR_LEN_MAX (9999)
#define HEADER "#{data#"
#define FOOTER "#}data#"
#define END_MSG "\n#file_completed#\n"

#define ASSERT assert
#define Fclose(x)  (x ? fclose(x),(x) = NULL : NULL)

