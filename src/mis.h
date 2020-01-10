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
 if(isFrame(X))
  setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));

 UNPROTECT(1);
 return(Ans);
}

enum nrm_mode {nmNone=0,nmSym=1,nmDirected=2};

SEXP miMatrix(SEXP X,SEXP Diag,SEXP Threads,enum nrm_mode mode){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,Threads,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
 SEXP Ans=PROTECT(allocMatrix(REALSXP,m,m));
 int zd=LOGICAL(Diag)[0];
 int *cXc=(int*)R_alloc(sizeof(int),2*n*nt);
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cA=cXc+(tn*n),*cB=cXc+((nt+tn)*n),da;
  struct ht *ht=hta[tn];
  #pragma omp for schedule(static,4)
  for(int a=0;a<m;a++){
   da=0;
   for(int b=0;b<=a;b++){
    if(a==b && zd){
     score[a*m+b]=0.;
     continue;
    }
    fillHt(ht,n,nx[a],x[a],nx[b],x[b],NULL,da?NULL:cA,cB,0);da=1;
    if(mode==nmNone)
     score[a*m+b]=score[b*m+a]=miHt(ht,cA,cB);
    else if(mode==nmSym)
     score[a*m+b]=score[b*m+a]=nmiHt(ht,cA,cB);
    else if(mode==nmDirected){
     double mi=miHt(ht,cA,cB);
     //TODO: Is this a right order?
     score[a*m+b]=mi/hC(n,nx[a],cA);
     score[b*m+a]=mi/hC(n,nx[b],cB);
    }
   }
  }
 }

 //Copy attribute names into both dimensions
 if(isFrame(X)){
  SEXP dimnames=PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(dimnames,0,getAttrib(X,R_NamesSymbol));
  SET_VECTOR_ELT(dimnames,1,getAttrib(X,R_NamesSymbol));
  setAttrib(Ans,R_DimNamesSymbol,dimnames);
  UNPROTECT(1);
 }

 UNPROTECT(1);
 return(Ans);
}

SEXP C_miMatrix(SEXP X,SEXP Diag,SEXP Threads){
 return(miMatrix(X,Diag,Threads,nmNone));
}

SEXP C_nmiMatrix(SEXP X,SEXP Diag,SEXP Threads){
 return(miMatrix(X,Diag,Threads,nmSym));
}

SEXP C_dnmiMatrix(SEXP X,SEXP Diag,SEXP Threads){
 return(miMatrix(X,Diag,Threads,nmDirected));
}
