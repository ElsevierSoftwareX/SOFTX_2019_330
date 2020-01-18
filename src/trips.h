int n2l(int n){
 if(n<2) return(0);
 //Simply n*(n-1)*(n-2)/6, 
 //but with an overflow avoidance
 int r=n%6;
 int q=(n-1)*(n-2);
 return((n/6)*q+(r*q)/6);
}

int ij2q(int i,int j){
 return(j*(j-1)/2+i);
}

SEXP C_tri(SEXP X,SEXP Threads){
 int n,m,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,R_NilValue,R_NilValue,Threads,&hta,&n,&m,NULL,NULL,NULL,&x,&nx,&nt);
 if(m<3) error("Cannot process less than three columns");
 if(m>2345) error("Too many features (>2345)");


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

 //Number of pairs for MI cache
 int np=m*(m-1)/2; //Won't overflow i32 since m<=2345
 double *mi_cache=(double*)R_alloc(sizeof(double),np);

 //Other caches
 int *caches=(int*)R_alloc(sizeof(int),nt*n*3);

 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
  int *cA  =caches+3*n*tn+0*n,
      *cB  =caches+3*n*tn+1*n,
      *ab  =caches+3*n*tn+2*n;
 
  //First step; fill the MI cache
  #pragma omp for schedule(static,4)
  for(int ei=0;ei<m;ei++){
   int da=0;
   for(int ej=ei+1;ej<m;ej++){
    fillHt(ht,n,nx[ei],x[ei],nx[ej],x[ej],NULL,da?NULL:cA,cB,0);da=1;
    mi_cache[ij2q(ei,ej)]=miHt(ht,cA,cB);
   }
  }

  int *cC=cA,*cAB=cB;
  int oa=-1,ob=-1,nab;
  //Second step; use I(A;B;C)=I(A;C)+I(B;C)-I(AB;C)
  // first two from MI cache, last by holding AB mix
  #pragma omp for
  for(int e=0;e<l;e++){
   int i=ii[e]-1,j=ji[e]-1,k=ki[e]-1;
   int *a=x[i],*b=x[j],*c=x[k];
   int na=nx[i],nb=nx[j],nc=nx[k];

   if(oa!=i || ob!=j){
    //New AB; we cache ab, nab and cAB
    nab=fillHt(ht,n,na,a,nb,b,ab,NULL,NULL,1);
    fillHt(ht,n,nc,c,nab,ab,NULL,cC,cAB,0);
    oa=i;ob=j;
   }else{
    fillHt(ht,n,nc,c,nab,ab,NULL,cC,NULL,0);
   }
   v[e]=mi_cache[ij2q(i,k)]+mi_cache[ij2q(j,k)]-miHt(ht,cC,cAB);
  }
 }

 UNPROTECT(6);
 return(Ans);
}
