SEXP jmi(SEXP X,SEXP Y,SEXP Z,SEXP Threads,int nrm){
 int n,m,ny,*y,nz,*z,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 if(length(Z)!=n) error("Z vector size mismatch");
 z=convertSEXP(*hta,n,Z,&nz);

 int *cXZc=(int*)R_alloc(sizeof(int),n*nt),
  *cY=(int*)R_alloc(sizeof(int),n),
  *xzc=(int*)R_alloc(sizeof(int),n*nt);

 for(int e=0;e<ny;e++) cY[e]=0;
 for(int e=0;e<n;e++) cY[y[e]-1]++;
 
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cXZ=cXZc+(tn*n),*xz=xzc+(tn*n);
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   //Mix X and Z
   int nxz=fillHt(ht,n,nz,z,nx[e],x[e],xz,NULL,NULL,1);
   fillHt(ht,n,ny,y,nxz,xz,NULL,NULL,cXZ,0);
   if(nrm) score[e]=nmiHt(ht,cY,cXZ); 
    else score[e]=miHt(ht,cY,cXZ);
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_jmi(SEXP X,SEXP Y,SEXP Z,SEXP Threads){
 return(jmi(X,Y,Z,Threads,0));
}

SEXP C_njmi(SEXP X,SEXP Y,SEXP Z,SEXP Threads){
 return(jmi(X,Y,Z,Threads,1));
}

SEXP C_cmi(SEXP X,SEXP Y,SEXP Z,SEXP Threads){
 int n,m,ny,*y,nz,*z,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 if(length(Z)!=n) error("Z vector size mismatch");
 z=convertSEXP(*hta,n,Z,&nz);
 int *cXZc=(int*)R_alloc(sizeof(int),n*nt);
 int *xzc=(int*)R_alloc(sizeof(int),n*nt);

 int *cZ=(int*)R_alloc(sizeof(int),n);
 int *yz2z=(int*)R_alloc(sizeof(int),n);
 int *yz=(int*)R_alloc(sizeof(int),n);
 int *cYZ=(int*)R_alloc(sizeof(int),n);
 //yz starts from one
 int nyz=fillHt(hta[0],n,ny,y,nz,z,yz,NULL,cZ,1);
 mixCountsHt(hta[0],cYZ);
 transHt(hta[0],NULL,yz2z);
 
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);
 for(int e=0;e<m;e++) score[e]=0.;

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  int *cXZ=cXZc+tn*n,*xz=xzc+tn*n;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   int nxz=fillHt(ht,n,nx[e],x[e],nz,z,xz,NULL,NULL,1);
   fillHt(ht,n,nxz,xz,nyz,yz,NULL,cXZ,NULL,0);
   score[e]=cmiHt(ht,cXZ,cYZ,yz2z,cZ);
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_cmiMatrix(SEXP X,SEXP Z,SEXP Diag,SEXP Threads){
 int n,m,*nx,**x,nt,nz,*z;
 struct ht **hta;
 int zd=LOGICAL(Diag)[0];
 prepareInput(X,Z,R_NilValue,Threads,&hta,&n,&m,NULL,&z,&nz,&x,&nx,&nt);
 SEXP Ans=PROTECT(allocMatrix(REALSXP,m,m));

 //Space for X_iW, essentially a second copy of X
 int *zx=(int*)R_alloc(sizeof(int),n*m);
 int *cZ=(int*)R_alloc(sizeof(int),n);
 for(int e=0;e<nz;e++) cZ[e]=0;
 for(int e=0;e<n;e++) cZ[z[e]-1]++;
 int *cXc=(int*)R_alloc(sizeof(int),2*n*nt);
 int *nzx=(int*)R_alloc(sizeof(int),m);
 //Collapsed state translation matrix
 int *xz2z=(int*)R_alloc(sizeof(int),n*m);
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   nzx[e]=fillHt(ht,n,nx[e],x[e],nz,z,zx+(n*e),NULL,NULL,1);
   transHt(ht,NULL,xz2z+n*e);
  }
  

  int *cAZ=cXc+(tn*n),*cBZ=cXc+((nt+tn)*n),da;
  #pragma omp for schedule(static,4)
  for(int a=0;a<m;a++){
   da=0;
   for(int b=0;b<=a;b++){
    if(a==b && zd){
     score[a*m+b]=0.;
     continue;
    }
    fillHt(ht,n,nzx[b],zx+(n*b),nzx[a],zx+(n*a),NULL,cBZ,da?NULL:cAZ,0);da=1;
    score[b*m+a]=score[a*m+b]=cmiHt(ht,cBZ,cAZ,xz2z+n*a,cZ);
   }
  }
 }

 //Copy attribute names into both dimensions
 SEXP dimnames=PROTECT(allocVector(VECSXP,2));
 SET_VECTOR_ELT(dimnames,0,getAttrib(X,R_NamesSymbol));
 SET_VECTOR_ELT(dimnames,1,getAttrib(X,R_NamesSymbol));
 setAttrib(Ans,R_DimNamesSymbol,dimnames);
 
 UNPROTECT(2);
 return(Ans);
}

SEXP jmiMatrix(SEXP X,SEXP W,SEXP Diag,SEXP Threads,int nrm){
 int zd=LOGICAL(Diag)[0];
 int n,m,*nx,**x,nt,nw,*w;
 struct ht **hta;
 prepareInput(X,W,R_NilValue,Threads,&hta,&n,&m,NULL,&w,&nw,&x,&nx,&nt);
 SEXP Ans=PROTECT(allocMatrix(REALSXP,m,m));

 //Space for X_iW, essentially a second copy of X
 int *wx=(int*)R_alloc(sizeof(int),n*m);

 int *cXc=(int*)R_alloc(sizeof(int),2*n*nt);
 int *nwx=(int*)R_alloc(sizeof(int),m);
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++)
   nwx[e]=fillHt(ht,n,nx[e],x[e],nw,w,wx+(n*e),NULL,NULL,1);
  

  //Calculate I(A,WB)
  int *cA=cXc+(tn*n),*cB=cXc+((nt+tn)*n),da;
  #pragma omp for schedule(static,4)
  for(int a=0;a<m;a++){
   da=0;
   for(int b=0;b<m;b++){
    if(a==b && zd){
     score[a*m+b]=0.;
     continue;
    }
    fillHt(ht,n,nx[a],x[a],nwx[b],wx+(n*b),NULL,da?NULL:cA,cB,0);da=1;
    if(nrm)
     score[a*m+b]=nmiHt(ht,cA,cB); else
     score[a*m+b]=miHt(ht,cA,cB);
   }
  }
 }

 //Copy attribute names into both dimensions
 SEXP dimnames=PROTECT(allocVector(VECSXP,2));
 SET_VECTOR_ELT(dimnames,0,getAttrib(X,R_NamesSymbol));
 SET_VECTOR_ELT(dimnames,1,getAttrib(X,R_NamesSymbol));
 setAttrib(Ans,R_DimNamesSymbol,dimnames);
 
 UNPROTECT(2);
 return(Ans);
}

SEXP C_jmiMatrix(SEXP X,SEXP W,SEXP Diag,SEXP Threads){
 return(jmiMatrix(X,W,Diag,Threads,0));
}

SEXP C_njmiMatrix(SEXP X,SEXP W,SEXP Diag,SEXP Threads){
 return(jmiMatrix(X,W,Diag,Threads,1));
}

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
