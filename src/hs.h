SEXP C_h(SEXP X,SEXP Threads){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,Threads,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  //Re-purpose ht memory
  int *cX=(int*)((hta[tn])->cnt);
  #pragma omp for
  for(int e=0;e<m;e++){
   for(int ee=0;ee<nx[e];ee++) cX[ee]=0;
   for(int ee=0;ee<n;ee++) cX[x[e][ee]-1]++;

   double H=0.;
   for(int ee=0;ee<nx[e];ee++) if(cX[ee]>0){
    double cAB=((double)(cX[ee]));
    H+=-cAB*log(cAB/((double)n));
   }
   score[e]=H/((double)n);
  }
 }
 //Copy attribute names
 if(isFrame(X))
  setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}

SEXP C_jh(SEXP X,SEXP Y,SEXP Threads){
 int n,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,R_NilValue,Threads,&hta,&n,&m,NULL,&y,&ny,&x,&nx,&nt);

 int *cache=(int*)R_alloc(sizeof(int),n*nt*2);

 SEXP Ans=PROTECT(allocVector(REALSXP,m));
 double *score=REAL(Ans);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num(),*cXY=cache+(tn*n*2),*xy=cache+(tn*n*2)+n;
  struct ht *ht=hta[tn];
  #pragma omp for
  for(int e=0;e<m;e++){
   //Mix X and Y
   int nxy=fillHt(ht,n,ny,y,nx[e],x[e],xy,NULL,NULL,0);

   //Count states
   for(int ee=0;ee<nxy;ee++) cXY[ee]=0;
   for(int ee=0;ee<n;ee++) cXY[xy[ee]]++;

   double H=0.;
   for(int ee=0;ee<nxy;ee++) if(cXY[ee]>0){
    double cAB=((double)(cXY[ee]));
    H+=-cAB*log(cAB/((double)n));
   }
   score[e]=H/((double)n);
  }
 }
 //Copy attribute names
 if(isFrame(X))
  setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));

 UNPROTECT(1);
 return(Ans);
}

