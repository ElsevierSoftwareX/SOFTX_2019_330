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
 setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 
 UNPROTECT(1);
 return(Ans);
}
