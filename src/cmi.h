SEXP C_CMI(SEXP X,SEXP Y,SEXP K,SEXP Threads){
 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 double bs=0.; int *cY,*ctmp,bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Save selected X as W and discard from further consideration
 int *w=(int*)R_alloc(sizeof(int),n);
 for(int e=0;e<n;e++) w[e]=x[bi][e];
 int nw=nx[bi]; x[bi]=NULL;

 int *w_tmp=(int*)R_alloc(sizeof(int),n);

 //Yet put it as a first W
 int *yw=(int*)R_alloc(sizeof(int),n);
 int *cW=(int*)R_alloc(sizeof(int),n);
 int *cYW=(int*)R_alloc(sizeof(int),n);
 int *yw2w=(int*)R_alloc(sizeof(int),n);
 int nyw=fillHt(hta[0],n,ny,y,nw,w,yw,NULL,cW,1);
 mixCountsHt(hta[0],cYW);
 transHt(hta[0],NULL,yw2w);

 //And as a first selected
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;

 //Time for an actual algorithm
 bs=0.;
 int *cXWc=(int*)R_alloc(sizeof(int),n*nt);
 int *xwc=(int*)R_alloc(sizeof(int),n*nt);

 int ke=k;
 #pragma omp parallel num_threads(nt)
 for(int e=1;e<ke;e++){
  double tbs=0.;
  int tbi=-1,tn=omp_get_thread_num();
  int *xw=xwc+n*tn,*cXW=cXWc+n*tn;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   int nxw=fillHt(ht,n,nx[ee],x[ee],nw,w,xw,NULL,NULL,1);
   fillHt(ht,n,nxw,xw,nyw,yw,NULL,cXW,NULL,0);
   double cs=cmiHt(ht,cXW,cYW,yw2w,cW);
   if(cs>tbs){
    tbs=cs; tbi=ee;
   }
  }
  #pragma omp critical
  if((tbs>bs) || (tbs==bs && tbi<bi)){
   bs=tbs;
   bi=tbi;
  }

  #pragma omp barrier
  #pragma omp single
  {
   if(bs>0){
    nw=fillHt(ht,n,nw,w,nx[bi],x[bi],w_tmp,NULL,NULL,1); w=w_tmp;
    mixCountsHt(ht,cW);
    nyw=fillHt(ht,n,ny,y,nw,w,yw,cY,cW,1);
    mixCountsHt(ht,cYW);
    transHt(ht,NULL,yw2w);
    x[bi]=NULL;

    score[e]=bs; idx[e]=bi+1; bs=0.;
   }else ke=e;
  }
 }

 Ans=finishAns(ke,Ans,X);

 UNPROTECT(1);
 return(Ans);
}
