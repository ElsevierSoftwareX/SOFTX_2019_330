SEXP C_mi(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);
 int *cXc=(int*)R_alloc(sizeof(int),n*nt);
 int *cYc=(int*)R_alloc(sizeof(int),n*nt);
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cX=cXc+(tn*n),*cY=cYc+(tn*n),dy=0;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   fillHt(ht,n,ny,y,nx[e],x[e],NULL,dy?NULL:cY,cX,0); dy=1;
   score[e]=miHt(ht,cY,cX);
  }
 }
 //Copy attribute names
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_miMatrix(SEXP X,SEXP Diag,SEXP Threads){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,Threads,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
 SEXP Ans=PROTECT(allocMatrix(REALSXP,m,m));
 int zd=LOGICAL(Diag)[0];
 int *cXc=(int*)R_alloc(sizeof(int),2*n*nt);
 double *score=REAL(Ans);

 //FIXME: Inefficient re-calculation of cA, also work spread over threads

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cA=cXc+(tn*n),*cB=cXc+((nt+tn)*n);
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m*m;e++){
   int a=e%m,b=e/m;
   if(a<b) continue;
   if(a==b && zd){
    score[e]=0.;
    continue;
   }
   fillHt(ht,n,nx[a],x[a],nx[b],x[b],NULL,cA,cB,0);
   score[a*m+b]=score[b*m+a]=miHt(ht,cA,cB);
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

