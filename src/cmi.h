enum cmi_jmi_mode {cjmCMI=791,cjmJMI=792,cjmNJMI=793};

//TODO: Drop CMI stuff from here, and migrate this to separate file, jmi.h
SEXP C_cmi_jmi(SEXP X,SEXP Y,SEXP Z,SEXP Mode,SEXP Threads){
 if(length(Mode)!=1) error("Invalid mode");
 int mode=INTEGER(Mode)[0];
 if(mode!=cjmCMI && mode!=cjmJMI && mode!=cjmNJMI)
  error("Invalid mode"); 

 int n,m,ny,*y,nz,*z,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 if(length(Z)!=n) error("Z vector size mismatch");
 z=convertSEXP(*hta,n,Z,&nz);

 int *cXZc=(int*)R_alloc(sizeof(int),n*nt),
  *cY=(int*)R_alloc(sizeof(int),n),
  *xzc=(int*)R_alloc(sizeof(int),n*nt);

 double scoreOff;
 {
  //Calculate I(Y;Z), which is a constant factor for CMI
  int *cZ=cXZc;
  //This is always needed for cY
  fillHt(*hta,n,ny,y,nz,z,NULL,cY,cZ,0);
  scoreOff=(mode==cjmCMI)?-miHt(*hta,cY,cZ):0.;
 }
 
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
   // I(X;Y|Z)=I(Y;X,Z)-I(Y;Z)
   if(mode!=cjmNJMI){
    //CMI or JMI=I(Y;X,Z)
    score[e]=miHt(ht,cY,cXZ)+scoreOff;
   }else{
    //NJMI (i.e. DISR-like score)=I(Y;X,Z)/H(X,Y,Z)
    score[e]=nmiHt(ht,cY,cXZ);
   }
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_cmi(SEXP X,SEXP Y,SEXP Z,SEXP Threads){
 int n,m,ny,*y,nz,*z,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 if(length(Z)!=n) error("Z vector size mismatch");
 z=convertSEXP(*hta,n,Z,&nz);
 int *cXYZc=(int*)R_alloc(sizeof(int),n*nt);
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
  int *cXYZ=cXYZc+tn*n,*cXZ=cXZc+tn*n,*xz=xzc+tn*n;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   int nxz=fillHt(ht,n,nx[e],x[e],nz,z,xz,NULL,NULL,1);
   fillHt(ht,n,nxz,xz,nyz,yz,NULL,cXZ,cXYZ,0);
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

SEXP C_jmiMatrix(SEXP X,SEXP W,SEXP Diag,SEXP Threads){
 int n,m,*nx,**x,nt,nw,*w;
 struct ht **hta;
 prepareInput(X,W,R_NilValue,Threads,&hta,&n,&m,NULL,&w,&nw,&x,&nx,&nt);
 SEXP Ans=PROTECT(allocMatrix(REALSXP,m,m));

 //Space for X_iW, essentially a second copy of X
 int *wx=(int*)R_alloc(sizeof(int),n*m);

 int zd=LOGICAL(Diag)[0];
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

