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

