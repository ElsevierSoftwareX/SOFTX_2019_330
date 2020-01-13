SEXP C_join(SEXP X){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,R_NilValue,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
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
 UNPROTECT(1);

 return(Ans);
}
