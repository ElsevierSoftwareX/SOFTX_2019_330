SEXP C_EJMI(SEXP X,SEXP Y,SEXP K,SEXP Threads){
 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 double bs=0.; int *cY,*ctmp,bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Save selected X as W and discard from further consideration
 int* w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1;
 
 //Time for an actual algorithm
 double *as=(double*)R_alloc(sizeof(double),m); //Accumulated score
 for(int e=0;e<m;e++) as[e]=0.;
 int *wxc=(int*)R_alloc(sizeof(int),n),*cWXc=ctmp,ke=k;
 bs=0.;

 for(int e=1;e<ke;e++){
  double tbs=0.;
  int tbi=-1,tn=0;
  struct ht *ht=hta[tn];
  int *wx=wxc,*cWX=cWXc;
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   //Mix x[ee] with lx making wx
   int nwx=fillHt(ht,n,nx[ee],x[ee],nw,w,wx,NULL,NULL,1);

   //Make MI of mix and Y and increase its accumulated score
   fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0); //cY stuff is red.
   as[ee]+=miHt(ht,cY,cWX);

   if(as[ee]>tbs){
    tbs=as[ee]; tbi=ee;
   }
  }
  if(tbs>bs){
   bs=tbs;
   bi=tbi;
  }
  
  if(bs>0.){
   w=x[bi]; nw=nx[bi]; x[bi]=NULL; 
   score[e]=bs; idx[e]=bi+1;
   bs=0.;
  }else ke=e;
 }

 Ans=finishAns(ke,Ans,X);

 UNPROTECT(1);
 return(Ans);
}

