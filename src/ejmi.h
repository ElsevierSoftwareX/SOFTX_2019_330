#include "pcg.h"

SEXP C_EJMI(SEXP X,SEXP Y,SEXP K,SEXP Iters,SEXP Prob,SEXP Threads){
 if(length(Iters)!=1) error("Iteration count should be a single value");
 int iters=INTEGER(Iters)[0];
 if(iters<1) error("Iteration count should be a positive integer");

 if(length(Prob)!=1) error("Acceptance probability should be a single value");
 double p=REAL(Prob)[0];
 if(p<=0. || p>=1.) error("Acceptance probability must be in (0;1)");
 uint32_t thresh=((double)(~(uint32_t)0))*p;
 uint64_t seed=seed_from_r();

 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 //Prepare a place for the answer
 SEXP Ans; PROTECT(Ans=allocVector(INTSXP,m));
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 int *ans=INTEGER(Ans);
 for(int e=0;e<m;e++) ans[e]=0;

 //Scores per feature; one vector per thread
 double *scores=(double*)R_alloc(sizeof(double),m*nt);
 //MI scores per feature; one for all
 double *mi=(double*)R_alloc(sizeof(double),m);
 //cY vector, common for all threads
 int *cY=(int*)R_alloc(sizeof(int),ny);
 for(int e=0;e<ny;e++) cY[e]=0;
 for(int e=0;e<n;e++) cY[y[e]-1]++;

 //Working place for each thread
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *wxc=(int*)R_alloc(sizeof(int),n*nt);
 int *masks=(int*)R_alloc(sizeof(int),m*nt);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  struct ht *ht=hta[tn];

  //First, we calculate MI(X,Y)
  int *cX=cXc+tn*n;
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,NULL,cX,0);
   mi[e]=miHt(ht,cY,cX);
  }
  //Implicit barrier here, mi is ready
  
  //Re-init 
  int *used=masks+tn*m,*cWX=cX,*wx=wxc+tn*n;
  double *score=scores+tn*m;
  for(int e=0;e<m;e++) used[e]=-1;

  //Loop over ensemble members
  #pragma omp for
  for(int em=0;em<iters;em++){
   //Clear the mask of used and feature score accumulator
   for(int e=0;e<m;e++) score[e]=0.;

   //Init the member-local RNG state;
   // waste one number for warm-up
   uint64_t rngs=seed*em; rng(&rngs,em);

   //Select first feature from the pre-computed mi vector
   double bs=0.; int bi=-1;
   for(int e=0;e<m;e++)
    if(mi[e]>=bs && rng(&rngs,em)<thresh){
     bs=mi[e]; bi=e;
    }

   if(bi==-1) continue; //No stem, no fun
   used[bi]=em;

   #pragma omp critical
   {
    ans[bi]++;
   }
   int *w=x[bi],nw=nx[bi];

   for(int e=1;e<k;e++){
    bs=0;
    for(int ee=0;ee<m;ee++){
     //Ignore already integrated features
     if(used[ee]==em) continue;

     //Mix x[ee] with lx making wx
     int nwx=fillHt(ht,n,nx[ee],x[ee],nw,w,wx,NULL,NULL,1);

     //Make MI of mix and Y and increase its accumulated score
     fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0);
     score[ee]+=miHt(ht,cY,cWX);

     //TODO: >= or >... Which is better?
     if(score[ee]>=bs && rng(&rngs,em)<thresh){
      bs=score[ee];
      bi=ee;
     }
    }

    if(bs>0.){
     //Give bi a hit
     #pragma omp critical
     {
      ans[bi]++;
     }
     w=x[bi]; nw=nx[bi]; used[bi]=em;
    }else break;
   }
  }
 }

 UNPROTECT(1);
 return(Ans);
}

SEXP C_EJMI2(SEXP X,SEXP Y,SEXP K,SEXP Iters,SEXP Prob,SEXP Threads){
 if(length(Iters)!=1) error("Iteration count should be a single value");
 int iters=INTEGER(Iters)[0];
 if(iters<1) error("Iteration count should be a positive integer");

 if(length(Prob)!=1) error("Acceptance probability should be a single value");
 double p=REAL(Prob)[0];
 if(p<=0. || p>=1.) error("Acceptance probability must be in (0;1)");
 uint32_t thresh=((double)(~(uint32_t)0))*p;
 uint64_t seed=seed_from_r();

 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 //Prepare a place for the answer
 SEXP Ans; PROTECT(Ans=allocVector(INTSXP,m));
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 int *ans=INTEGER(Ans);
 for(int e=0;e<m;e++) ans[e]=0;

 //Scores per feature; one vector per thread
 double *scores=(double*)R_alloc(sizeof(double),m*nt);
 //MI scores per feature; one for all
 double *mi=(double*)R_alloc(sizeof(double),m);
 //cY vector, common for all threads
 int *cY=(int*)R_alloc(sizeof(int),ny);
 for(int e=0;e<ny;e++) cY[e]=0;
 for(int e=0;e<n;e++) cY[y[e]-1]++;

 //Working place for each thread
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *wxc=(int*)R_alloc(sizeof(int),n*nt);
 int *masks=(int*)R_alloc(sizeof(int),m*nt);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  struct ht *ht=hta[tn];

  //First, we calculate MI(X,Y)
  int *cX=cXc+tn*n;
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,NULL,cX,0);
   mi[e]=miHt(ht,cY,cX);
  }
  //Implicit barrier here, mi is ready
  
  //Re-init 
  int *used=masks+tn*m,*cWX=cX,*wx=wxc+tn*n;
  double *score=scores+tn*m;
  for(int e=0;e<m;e++) used[e]=-1;

  //Loop over ensemble members
  #pragma omp for
  for(int em=0;em<iters;em++){
   //Clear the mask of used and feature score accumulator
   for(int e=0;e<m;e++) score[e]=0.;

   //Init the member-local RNG state;
   // waste one number for warm-up
   uint64_t rngs=seed*em; rng(&rngs,em);

   //Remove each feature with 1-p chance
   for(int e=0;e<m;e++)
    if(rng(&rngs,em)>thresh)
     used[e]=em;

   //Select first feature from the pre-computed mi vector
   double bs=0.; int bi=-1;
   for(int e=0;e<m;e++)
    if(used[e]!=em && mi[e]>=bs){
     bs=mi[e]; bi=e;
    }

   if(bi==-1) continue; //No stem, no fun
   used[bi]=em;

   #pragma omp critical
   {
    ans[bi]++;
   }
   int *w=x[bi],nw=nx[bi];

   for(int e=1;e<k;e++){
    bs=0;
    for(int ee=0;ee<m;ee++){
     //Ignore already integrated features
     if(used[ee]==em) continue;

     //Mix x[ee] with lx making wx
     int nwx=fillHt(ht,n,nx[ee],x[ee],nw,w,wx,NULL,NULL,1);

     //Make MI of mix and Y and increase its accumulated score
     fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0);
     score[ee]+=miHt(ht,cY,cWX);

     if(score[ee]>bs){
      bs=score[ee];
      bi=ee;
     }
    }

    if(bs>0.){
     //Give bi a hit
     #pragma omp critical
     {
      ans[bi]++;
     }
     w=x[bi]; nw=nx[bi]; used[bi]=em;
    }else break;
   }
  }
 }

 UNPROTECT(1);
 return(Ans);
}
