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
