
## compiling bam:

. compiling bam with an external lib, one can set the followings into MyConfig:

INCS += -I/home/alireza/Workstation/Elliptica_ID_Reader/include/ 
SPECIALLIBS += -L/home/alireza/Workstation/Elliptica_ID_Reader/lib/ -lelliptica_id_reader

## NOTEs and TODOs:
. the latest commit has not been tested for SBH and NSNS sections.

. for now, one can use the version v1.1 for SBH.

