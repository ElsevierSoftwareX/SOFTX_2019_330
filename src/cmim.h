SEXP C_CMIM(SEXP X,SEXP Y,SEXP K,SEXP Threads){
 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 double bs=0.,*ms=(double*)R_alloc(sizeof(double),m);
 int bi=0,*cYc,*ctmp;
 bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cYc,&ctmp,ms,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));
 
 int *lk=(int*)R_alloc(sizeof(int),m),
  **w=(int**)R_alloc(sizeof(int*),k),
  *nw=(int*)R_alloc(sizeof(int),k),
  *cW=(int*)R_alloc(sizeof(int),n*k),
  *wy=(int*)R_alloc(sizeof(int),n*k),
  *nwy=(int*)R_alloc(sizeof(int),k),
  *cWY=(int*)R_alloc(sizeof(int),n*k),
  *cXWc=(int*)R_alloc(sizeof(int),n*nt),
  *xwc=(int*)R_alloc(sizeof(int),n*nt),
  *wy2w=(int*)R_alloc(sizeof(int),n*k);
 w[0]=x[bi]; nw[0]=nx[bi]; x[bi]=NULL;
 for(int e=0;e<m;e++) lk[e]=0;
 double *score; int *idx,ke=k;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;

 //Commented out according to Fleuret's paper
 //for(int e=0;e<n;e++) ms[e]=INFINITY;
 
 #pragma omp parallel num_threads(nt)
 for(int e=1;e<ke;e++){
  #pragma omp single
  {
   bs=-INFINITY;
   int off=n*(e-1);
   nwy[e-1]=fillHt(hta[0],n,nw[e-1],w[e-1],ny,y,wy+off,cW+off,NULL,1);
   transHt(hta[0],wy2w+off,NULL);
   mixCountsHt(hta[0],cWY+off);
   //for(int ee=0;ee<nwy[e-1];ee++) (cWY+off)[ee]=hta[0]->cnt[ee].c;
  }
  double tbs=-INFINITY;
  int tbi=-1,tn=omp_get_thread_num(),*cX=cYc+(tn*n);
  int *xw=xwc+tn*n,*cXW=cXWc+tn*n;
  struct ht *ht=hta[tn];
  #pragma omp for schedule(dynamic)
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected or with no chance to be max
   if(!x[ee] || tbs>ms[ee]) continue;
   
   //Push forward towards e
   for(;lk[ee]<e;lk[ee]++){
    int ew=lk[ee],off=ew*n;
    int nxw=fillHt(ht,n,nx[ee],x[ee],nw[ew],w[ew],xw,cX,NULL,1);
    fillHt(ht,n,nxw,xw,nwy[ew],wy+off,NULL,cXW,NULL,0);
    double cmi=cmiHt(ht,cXW,cWY+off,wy2w+off,cW+off);

    ms[ee]=(ms[ee]<cmi)?ms[ee]:cmi;
    //Maybe it is already eliminated?
    if(tbs>ms[ee]) break;
   }

   //Check again
   if(ms[ee]>tbs){
    tbs=ms[ee]; tbi=ee;
   }
  }
  #pragma omp critical
  if((tbs>bs) || (tbs==bs && tbi<bi)){
   bs=tbs;
   bi=tbi;
  }
  #pragma omp barrier
  #pragma omp single
  if(bs>0.){
   w[e]=x[bi]; nw[e]=nx[bi]; x[bi]=NULL;
   score[e]=bs; idx[e]=bi+1;
  }else ke=e;
 }

 Ans=finishAns(ke,Ans,X);
 UNPROTECT(1);
 return(Ans);
}

