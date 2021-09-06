//TODO: Rename this file

SEXP C_max_jmi(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);


 int *cXZc=(int*)R_alloc(sizeof(int),n*nt),
  *cY=(int*)R_alloc(sizeof(int),n),
  *xzc=(int*)R_alloc(sizeof(int),n*nt);

 for(int e=0;e<ny;e++) cY[e]=0;
 for(int e=0;e<n;e++) cY[y[e]-1]++;

 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);
 for(int e=0;e<m;e++) score[e]=0.; //Required for bumping into max

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cXZ=cXZc+(tn*n),*xz=xzc+(tn*n);
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int a=1;a<m;a++){
   for(int b=0;b<a;b++){
    int nxz=fillHt(ht,n,nx[a],x[a],nx[b],x[b],xz,NULL,NULL,1);
    fillHt(ht,n,ny,y,nxz,xz,NULL,NULL,cXZ,0);
    double s=miHt(ht,cY,cXZ);
    if(s>score[a]) score[a]=s;
    if(s>score[b]) score[b]=s;
   }
  }
 }
 //Copy attribute names
 if(isFrame(X))
  setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));

 UNPROTECT(1);
 return(Ans);
}

//TODO: Adjust function name
SEXP C_max_cmi(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 int *cXZc=(int*)R_alloc(sizeof(int),n*nt),
  *cYZc=(int*)R_alloc(sizeof(int),n*nt),
  *cY=(int*)R_alloc(sizeof(int),n),
  *cZc=(int*)R_alloc(sizeof(int),n*nt),
  *yzc=(int*)R_alloc(sizeof(int),n*nt),
  *xzc=(int*)R_alloc(sizeof(int),n*nt),
  *yz2zc=(int*)R_alloc(sizeof(int),n*nt);


 for(int e=0;e<ny;e++) cY[e]=0;
 for(int e=0;e<n;e++) cY[y[e]-1]++;

 SEXP Ans=PROTECT(allocMatrix(REALSXP,2,m));
 double *score=REAL(Ans);
 for(int e=0;e<m;e++){
  score[2*e+1]=0.; //Required for bumping into max...
  score[2*e]=INFINITY; //...and into min
 }

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cXZ=cXZc+(tn*n),*xz=xzc+(tn*n),
      *cYZ=cYZc+(tn*n),*yz=yzc+(tn*n),*yz2z=yz2zc+(tn*n),*cZ=cZc+(tn*n);
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int ez=0;ez<m;ez++){
   int nyz=fillHt(ht,n,ny,y,nx[ez],x[ez],yz,NULL,cZ,1);
   mixCountsHt(ht,cYZ);
   transHt(ht,NULL,yz2z);
   for(int e=0;e<m;e++) if(ez!=e) {
    int nxz=fillHt(ht,n,nx[e],x[e],nx[ez],x[ez],xz,NULL,NULL,1);
    fillHt(ht,n,nxz,xz,nyz,yz,NULL,cXZ,NULL,0);
    double s=cmiHt(ht,cXZ,cYZ,yz2z,cZ);
    if(s>score[2*e+1]) score[2*e+1]=s;
    if(s<score[2*e]) score[2*e]=s; 
   }
  }
 }
 //Copy attribute names
 if(isFrame(X)){
  SEXP dimnames=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(dimnames,0,R_NilValue);
  SET_VECTOR_ELT(dimnames,1,getAttrib(X,R_NamesSymbol));
  setAttrib(Ans,R_DimNamesSymbol,dimnames);
  UNPROTECT(1);
 }

 UNPROTECT(1);
 return(Ans);
}
