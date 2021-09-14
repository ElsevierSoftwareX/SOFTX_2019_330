#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_thread_num() 0
 #define omp_get_max_threads() 1
 #define omp_set_num_threads(x)
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h> 

//Hash table

#include "ht.h"

//Shared functions

#include "shared.h"

//Feature selection algorithms

#include "cmi.h"
#include "cmim.h"
#include "mim.h"
#include "mrmr.h"
#include "disr.h"
#include "jmi.h"
#include "jmim.h"
#include "njmim.h"

//Ensemble selection

#include "ejmi.h"

//Gini impurity based

#include "jim.h"

//Feature scoring algorithms

#include "mis.h"
#include "cmis.h"
#include "jmis.h"
#include "maxjmis.h"
#include "ims.h"
#include "hs.h"
#include "trips.h"

//Kendall transformation

#include "kt.h"

//Auxiliary

#include "side.h"
#include "join.h"

//Registration

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_CallMethodDef R_CallDef[]={
 CALLDEF(C_engineTest,2),
 CALLDEF(C_getMi,2),
 CALLDEF(C_CMIM,4),
 CALLDEF(C_CMI,4),
 CALLDEF(C_JMI,4),
 CALLDEF(C_DISR,4),
 CALLDEF(C_JMIM,4),
 CALLDEF(C_NJMIM,4),
 CALLDEF(C_MIM,4),
 CALLDEF(C_MRMR,5),
 CALLDEF(C_mi,3),
 CALLDEF(C_miMatrix,3),
 CALLDEF(C_jmi,4),
 CALLDEF(C_njmi,4),
 CALLDEF(C_max_jmi,3),
 CALLDEF(C_max_cmi,3),
 CALLDEF(C_cmi,4),
 CALLDEF(C_cmiMatrix,4),
 CALLDEF(C_cmiMatrix2,3),
 CALLDEF(C_jmiMatrix,4),
 CALLDEF(C_njmiMatrix,4),
 CALLDEF(C_nmiMatrix,3),
 CALLDEF(C_dnmiMatrix,3),
 CALLDEF(C_im,3),
 CALLDEF(C_h,2),
 CALLDEF(C_jh,3),
 CALLDEF(C_JIM,4),
 CALLDEF(C_kt,1),
 CALLDEF(C_rkt,1),
 CALLDEF(C_join,1),
 CALLDEF(C_tri,2),
 CALLDEF(C_EJMI,6),
 CALLDEF(C_EJMI2,6),
 {NULL,NULL,0}
};

void attribute_visible R_init_praznik(DllInfo *dll){
 R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
 R_useDynamicSymbols(dll,FALSE);
 R_forceSymbols(dll,TRUE);
}
