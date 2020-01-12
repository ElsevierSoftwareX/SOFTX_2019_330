int n2l(int n){
 if(n<2) return(0);
 //Version avoiding overflow, quite nasty
// int ans=n;
// if(n%2){
//  ans*=(n-1)/2;
// }else{
//  ans=(ans/2)*(n-1);
// }
// if(ans%3){
//  ans*=(n-2)/3;
// }else{
//  ans=(ans/3)*(n-2);
// }
 return(n*(n-1)*(n-2)/6);
}

SEXP C_tri(SEXP X,SEXP Threads){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,Threads,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
 if(m<3) error("Cannot process less than three columns");

 int nl=n2l(m);
 SEXP I=PROTECT(allocVector(INTSXP,nl));
 SEXP J=PROTECT(allocVector(INTSXP,nl));
 SEXP K=PROTECT(allocVector(INTSXP,nl));
 SEXP V=PROTECT(allocVector(REALSXP,nl));
 int *ii=INTEGER(I),*ji=INTEGER(J),*ki=INTEGER(K);
 double *v=REAL(V);

 int l=0;
 for(int ei=0;ei<m;ei++)
  for(int ej=ei+1;ej<m;ej++)
   for(int ek=ej+1;ek<m;ek++){
    ii[l]=ei+1; ji[l]=ej+1; ki[l]=ek+1;
    l++;
   }

 SEXP Ans=PROTECT(allocVector(VECSXP,4));
 SET_VECTOR_ELT(Ans,0,I);
 SET_VECTOR_ELT(Ans,1,J);
 SET_VECTOR_ELT(Ans,2,K);
 SET_VECTOR_ELT(Ans,3,V);

 SEXP AnsN=PROTECT(NEW_CHARACTER(4));
 SET_STRING_ELT(AnsN,0,mkChar("Var1"));
 SET_STRING_ELT(AnsN,1,mkChar("Var2"));
 SET_STRING_ELT(AnsN,2,mkChar("Var3"));
 SET_STRING_ELT(AnsN,3,mkChar("MI"));
 setAttrib(Ans,R_NamesSymbol,AnsN);

 SEXP Names=getAttrib(X,R_NamesSymbol);
 setAttrib(I,R_LevelsSymbol,Names);
 setAttrib(I,R_ClassSymbol,mkString("factor"));
 setAttrib(J,R_LevelsSymbol,Names);
 setAttrib(J,R_ClassSymbol,mkString("factor"));
 setAttrib(K,R_LevelsSymbol,Names);
 setAttrib(K,R_ClassSymbol,mkString("factor"));

 int *caches=(int*)R_alloc(sizeof(int),nt*n*8);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
    int *cA  =caches+8*n*tn+0*n,
        *cB  =caches+8*n*tn+1*n,
        *cC  =caches+8*n*tn+2*n,
        *cAC =caches+8*n*tn+3*n,
        *cBC =caches+8*n*tn+4*n,
        *bc2c=caches+8*n*tn+5*n,
        *ac  =caches+8*n*tn+6*n,
        *bc  =caches+8*n*tn+7*n;
  #pragma omp for
  for(int e=0;e<l;e++){
   int i=ii[e]-1,j=ji[e]-1,k=ki[e]-1;
   int *a=x[i],*b=x[j],*c=x[k];
   int na=nx[i],nb=nx[j],nc=nx[k];

   //TODO: Cache something here, lots of things are 
   // calculated over and over.

   //int nab=
   fillHt(ht,n,na,a,nb,b,NULL,cA,cB,0);
   double mi=miHt(ht,cA,cB);

   int nbc=fillHt(ht,n,nb,b,nc,c,bc,NULL,cC,1);
   mixCountsHt(ht,cBC);
   transHt(ht,NULL,bc2c);

   int nac=fillHt(ht,n,na,a,nc,c,ac,NULL,NULL,1);
   fillHt(ht,n,nac,ac,nbc,bc,NULL,cAC,NULL,0);
   double cmi=cmiHt(ht,cAC,cBC,bc2c,cC);

   v[e]=mi-cmi;
  }
 }

 UNPROTECT(6);
 return(Ans);
}
