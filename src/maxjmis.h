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
   //#pragma omp for schedule(static,4)
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

