void vec2kt(SEXP V,int n,int *x){
 if(length(V)!=n) error("Invalid length");
 if(isFactor(V))
  if(!isOrdered(V) && length(getAttrib(V,R_LevelsSymbol))>2) error("Unordered factor with more than two levels");
 
 if(isInteger(V)||isLogical(V)||isFactor(V)||isOrdered(V)){ //Logical NA is NA_INTEGER
  int *v=INTEGER(V);
  for(int i=0;i<n;i++){
   if(v[i]==NA_INTEGER){
    for(int j=0;j<n;j++)
     if(i!=j) (x++)[0]=NA_INTEGER;
   }else{
    for(int j=0;j<n;j++) if(i!=j){
     if(v[j]==NA_INTEGER){
      (x++)[0]=NA_INTEGER;
     }else{
      (x++)[0]=(v[i]<=v[j])+2*(v[i]>=v[j]);
     }  
    }
   }
  }
  return;
 }

 if(isReal(V)){
  double *v=REAL(V);
  for(int i=0;i<n;i++){
   if(ISNAN(v[i])||v[i]==NA_REAL){
    for(int j=0;j<n;j++)
     if(i!=j) (x++)[0]=NA_INTEGER;
   }else{
    for(int j=0;j<n;j++) if(i!=j){
      if(ISNAN(v[j])||v[j]==NA_REAL){
       (x++)[0]=NA_INTEGER;
      }else{
       (x++)[0]=(v[i]<=v[j])+2*(v[i]>=v[j]);
      }  
     }
   }
  }
  return;
 }

 error("Unsupported input");
}

SEXP C_kt(SEXP X){
 SEXP Ktl; PROTECT(Ktl=NEW_CHARACTER(3));
 SET_STRING_ELT(Ktl,0,mkChar("<"));
 SET_STRING_ELT(Ktl,1,mkChar(">"));
 SET_STRING_ELT(Ktl,2,mkChar("="));

 if(isLogical(X)||isInteger(X)||isFactor(X)||isReal(X)){
  int n=length(X);
  if(n<2) error("Vector shorter than 2 elements");
  SEXP Ans; PROTECT(Ans=allocVector(INTSXP,n*(n-1)));
  vec2kt(X,n,INTEGER(Ans));
  setAttrib(Ans,R_LevelsSymbol,Ktl);
  setAttrib(Ans,R_ClassSymbol,mkString("factor"));
  UNPROTECT(2); //KtV + Levels
  return(Ans);
 }

 //TODO: Add lists one day?
 if(isFrame(X)){
  if(length(X)<1) error("Data frame with no columns");
  int n=length(VECTOR_ELT(X,0));
  if(n<2) error("Data frame with less than two rows");
  SEXP Ans; PROTECT(Ans=allocVector(VECSXP,length(X)));
  for(int e=0;e<length(X);e++){
   SEXP V=VECTOR_ELT(X,e);
   SEXP KtV; PROTECT(KtV=allocVector(INTSXP,n*(n-1)));
   vec2kt(V,n,INTEGER(KtV));
   setAttrib(KtV,R_LevelsSymbol,Ktl);
   setAttrib(KtV,R_ClassSymbol,mkString("factor"));
   SET_VECTOR_ELT(Ans,e,KtV);
   UNPROTECT(1); //KtV
  }
  //Copy names
  setAttrib(Ans,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
  UNPROTECT(2); //Ans and Ktl
  //Ans is a list, we want to use data.frame to make row names
  return(Ans);
 }
 error("Unsupported input");
}
