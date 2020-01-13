SEXP C_join(SEXP X){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,NULL,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
 struct ht *ht=hta[0];

 SEXP Ans=PROTECT(allocVector(INTSXP,n));
 int *z=INTEGER(Ans),*z2=(int*)R_alloc(sizeof(int),n);
 if(m%2) ISWAP(z,z2);
 int nz=1,nz2;
 for(int e=0;e<n;e++) z[e]=1;

 for(int e=0;e<m;e++){
  nz2=fillHt(ht,n,nz,z,nx[e],x[e],z2,NULL,NULL,1);
  ISWAP(z,z2);
  nz=nz2;
 }

 char buf[64];
 SEXP Levels=PROTECT(allocVector(STRSXP,nz));
 for(int e=0;e<nz;e++){
  snprintf(buf,64,"l%d",e+1);
  SET_STRING_ELT(Levels,e,mkChar(buf));
 }
 setAttrib(Ans,R_LevelsSymbol,Levels);
 setAttrib(Ans,R_ClassSymbol,mkString("factor"));


 UNPROTECT(2);
 return(Ans);
}
