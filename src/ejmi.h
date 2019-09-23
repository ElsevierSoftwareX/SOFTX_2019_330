#include "pcg.h"

SEXP C_EJMI(SEXP X,SEXP Y,SEXP K,SEXP Iters,SEXP Prob,SEXP Threads){
 //FIXME: Handle errors here
 int iters=INTEGER(Iters)[0];
 double p=REAL(Prob)[0];
 if(p<=0 || p>=1.) error("p must be in (0;1)");
 uint32_t thresh=((double)0xFFFFFFFF)*p;
 uint32_t seed=17; //FIXME: Get this from R's PRNG

 uint64_t rngs=seed;

 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 double bs=0.; int *cY,*ctmp,bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Save selected X as W and discard from further consideration
 int* w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 SEXP Ans; PROTECT(Ans=allocVector(INTSXP,m));
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 int *ans=INTEGER(Ans);
 for(int e=0;e<m;e++)
  ans[e]=0;
 
 //Time for an actual algorithm
 double *as=(double*)R_alloc(sizeof(double),m); //Accumulated score
 for(int e=0;e<m;e++) as[e]=0.;
 int *wx=(int*)R_alloc(sizeof(int),n),*cWX=ctmp;
 bs=0.;

 for(int e=1;e<k;e++){
  struct ht *ht=hta[0];
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   //Mix x[ee] with lx making wx
   int nwx=fillHt(ht,n,nx[ee],x[ee],nw,w,wx,NULL,NULL,1);

   //Make MI of mix and Y and increase its accumulated score
   fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0);
   as[ee]+=miHt(ht,cY,cWX);

   if(as[ee]>bs){
    uint32_t rrr=rng(&rngs,7);
    if(rrr<thresh){
     bs=as[ee]; bi=ee;
    }
   }
  }
  
  if(bs>0.){
   //Just give this arg a hit
   ans[bi]++;
   w=x[bi]; nw=nx[bi]; x[bi]=NULL; 
   bs=0.;
  }else continue;
 }

 UNPROTECT(1);
 return(Ans);
}

