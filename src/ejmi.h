#include "pcg.h"

SEXP C_EJMI(SEXP X,SEXP Y,SEXP K,SEXP Iters,SEXP Prob,SEXP Threads){
 //FIXME: Handle errors here
 int iters=INTEGER(Iters)[0];
 double p=REAL(Prob)[0];
 //if(p<=0 || p>=1.) error("p must be in (0;1)");
 if(p<0 || p>1.) error("p must be in (0;1)"); //FIXME: allowed for testing
 uint32_t thresh=((double)0xFFFFFFFF)*p;

//  int zn=1/p*100000,z9=0,z123=0;
//  uint64_t tr=17;
//  for(int e=0;e<zn;e++){
//   z9+=rng(&tr,9)<thresh;
//   z123+=rng(&tr,123)<thresh;
//  }
//  printf("peff=%0.3g %0.3g, p=%0.3g\n",
//   ((double)z9)/((double)zn),
//   ((double)z123)/((double)zn),
//   p
//  );

 //GetRNGState();
 uint32_t seed=17;//R_unif_index((double)0xFFFFFFFF); 
 //PutRNGState();

 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 //Prepare a place for the answer
 SEXP Ans; PROTECT(Ans=allocVector(INTSXP,m));
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 int *ans=INTEGER(Ans);
 for(int e=0;e<m;e++) ans[e]=0;

 //Scores per feature; one vector per thread plus one for mi
 double *scores=(double*)R_alloc(sizeof(double),m*nt+m);
 //Integer vectors, three per thread, for cY, cX/cWX and wx
 int *ints=(int*)R_alloc(sizeof(int),n*nt*3);
 //Mask vectors, to mark already used features
 int *masks=(int*)R_alloc(sizeof(int),m*nt);


 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),dY=0;
  struct ht *ht=hta[tn];

  //First, we calculate MI(X,Y)
  double *mi=scores+nt*m;
  int *cY=ints+tn*3*n,*cX=ints+(tn*3+1)*n;
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,dY?NULL:cY,cX,0); dY=1;
   mi[e]=miHt(ht,cY,cX);
  }
  //Implicit barrier here, mi and cY are ready
  
  //Re-init 
  int *used=masks+tn*m,*cWX=cX,*wx=ints+(tn*3+2)*n;
  double *score=scores+tn*m;
  for(int e=0;e<m;e++) used[e]=-1;

  //Loop over ensemble members
  #pragma omp for
  for(int em=0;em<iters;em++){
   //Clear the mask of used and feature score accumulator
   for(int e=0;e<m;e++) score[e]=0.;

   //Init the member-local RNG state
   uint64_t rngs=seed+em; rng(&rngs,em);

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

