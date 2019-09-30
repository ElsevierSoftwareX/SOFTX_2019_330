#include <stdint.h>

#define GET_A(x) ((x)>>32)
#define GET_B(x) ((x)&0x00000000ffffffff)
#define MAKE_AB(a,b) (((uint64_t)(a)<<32)|(b))

struct ht{
 struct hte **map;
 struct hte *cnt;
 int N;
 uint32_t nAB;
};

struct hte{
 uint64_t ab;
 struct hte *nxt;
 int c;
};

struct ht* R_allocHt(int N){
 struct ht *ans=(struct ht*)R_alloc(sizeof(struct ht),1);
 ans->N=N;
 ans->map=(struct hte**)R_alloc(sizeof(struct hte*),N);
 ans->cnt=(struct hte*)R_alloc(sizeof(struct hte),N);
 return(ans);
}

uint32_t static inline fillHt(struct ht *Q,int N,int nA,int *a,int nB,int *b,int *mix,int *cA,int *cB,int mixOff){
 if(cA) for(int e=0;e<nA;e++) cA[e]=0;
 if(cB) for(int e=0;e<nB;e++) cB[e]=0;
 int Neff=nA*nB; 

 if(Neff<N && !mix){
  //HT not needed, we use it as lookup table 
  for(int e=0;e<Neff;e++) Q->cnt[e].c=0;
  for(int e=0;e<N;e++){
   uint32_t _a=a[e]-1,_b=b[e]-1,hab=_a+_b*nA;
   uint64_t _ab=MAKE_AB(_a,_b);
   Q->cnt[hab].c++;
   Q->cnt[hab].ab=_ab;
   if(cA) cA[_a]++;
   if(cB) cB[_b]++;
  }
  return(Q->nAB=Neff);
 }

 if(Neff<N && mix){
  //LT, but map used to collapse mix
  uint32_t nAB=0;
  for(int e=0;e<Neff;e++) Q->map[e]=NULL;
  for(int e=0;e<N;e++){
   uint32_t _a=a[e]-1,_b=b[e]-1,hab=_a+_b*nA;
   uint64_t _ab=MAKE_AB(_a,_b);
   struct hte **E=Q->map+hab;
   if(!*E){
    Q->cnt[nAB].ab=_ab;
    Q->cnt[nAB].c=1;
    *E=Q->cnt+nAB;
    nAB++;
   }else (*E)->c++;

   if(cA) cA[_a]++;
   if(cB) cB[_b]++;
   mix[e]=(*E-Q->cnt)+mixOff;
  }
  return(Q->nAB=nAB);
 }
 
 //Regular HT
 uint32_t nAB=0;
 for(int e=0;e<N;e++) Q->map[e]=NULL;
 
 for(int e=0;e<N;e++){
  uint32_t _a=a[e]-1,_b=b[e]-1,hab=(_a^_b)%N;
  uint64_t _ab=MAKE_AB(_a,_b);

  struct hte **E;
  for(E=Q->map+hab;(*E)&&(*E)->ab!=_ab;E=&((*E)->nxt));

  if(!*E){
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   Q->cnt[nAB].c=1;
   *E=Q->cnt+nAB;
   nAB++;
  }else (*E)->c++;

  if(cA) cA[_a]++;
  if(cB) cB[_b]++;
  if(mix) mix[e]=(*E-Q->cnt)+mixOff;
 }
 return(Q->nAB=nAB);
}

uint32_t static inline fillHtOne(struct ht *Q,int N,int *in,int *out,int mixOff){
 uint32_t nAB=0;
 for(int e=0;e<N;e++) Q->map[e]=NULL;
 for(int e=0;e<N;e++){
  if(in[e]==NA_INTEGER) error("NA values are not allowed");
  uint64_t _ab=(uint64_t)(in[e]);
  uint32_t hab=_ab%N;

  struct hte **E;
  for(E=Q->map+hab;E[0] && E[0]->ab!=_ab;E=&(E[0]->nxt));

  if(!E[0]){
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   E[0]=Q->cnt+nAB;
   nAB++;
  }
  out[e]=(E[0]-Q->cnt)+mixOff;
 }
 return(nAB);
}

double miHt(struct ht *Q,int *cA,int *cB){
 double ans=0.,N=Q->N;
 for(int e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double cAB=Q->cnt[e].c,
   _cA=cA[GET_A(Q->cnt[e].ab)],
   _cB=cB[GET_B(Q->cnt[e].ab)];
  ans+=cAB*log((cAB*N)/(_cA*_cB));
 }
 return(ans/N);
}

double nmiHt(struct ht *Q,int *cA,int *cB){
 double I=0.,H=0.,N=Q->N;
 for(int e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double cAB=Q->cnt[e].c,
   _cA=cA[GET_A(Q->cnt[e].ab)],
   _cB=cB[GET_B(Q->cnt[e].ab)];
  I+=cAB*log((cAB*N)/(_cA*_cB));
  H+=-cAB*log(cAB/N);
 }
 // NaN when I=H=0 is expected
 return(I/H);
}

//Note translation vectors are zero-based, i.e. one must subtract
// mixOff like mix2a[mix[e]-mixOff]-> a (in original offset)
void transHt(struct ht *Q,int *mix2a,int *mix2b){
 for(int e=0;e<Q->nAB;e++){
  if(mix2a)
   mix2a[e]=GET_A(Q->cnt[e].ab);
  if(mix2b)
   mix2b[e]=GET_B(Q->cnt[e].ab);
 }
}

void mixCountsHt(struct ht *Q,int *c){
 for(int e=0;e<Q->nAB;e++)
  c[e]=Q->cnt[e].c;
}

double cmiHt(struct ht *Q,int *cXW,int *cYW,int *yw2w,int *cW){
 double I=0.,N=Q->N;
 for(int e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double _cXYW=((double)Q->cnt[e].c),
   _cXW=cXW[GET_A(Q->cnt[e].ab)],
   _cYW=cYW[GET_B(Q->cnt[e].ab)],
   _cW=cW[yw2w[GET_B(Q->cnt[e].ab)]];
  I+=_cXYW*log(_cXYW*_cW/_cXW/_cYW);
 }
 return(I/N);
}

double hC(int n,int nc,int *c){
 double H=0.,N=n;
 for(int e=0;e<nc;e++) 
  if(c[e]){
   double C=c[e];
   H+=-C*log(C/N);
  }
 return(H/((double)n));
}

//Impurity is calculated in two parts,
// p_xy^2/p_x^2 ~ x and constant p_y^2 ...
double imHt(struct ht *Q,int *cA){
 double ans=0.,N=Q->N;
 for(int e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double cAB=Q->cnt[e].c,
   _cA=cA[GET_A(Q->cnt[e].ab)];
  ans+=cAB*cAB/_cA;
 }
 return(ans/N);
}
//...here
double imOff(int *cB,int nB,int N){
 double ans=0.;
 for(int e=0;e<nB;e++){
  double pa=((double)cB[e])/((double)N);
  ans+=pa*pa;
 }
 return(ans);
}
//...or here, if one only has a raw y vector
double imOffRaw(struct ht *Q,int N,int *b,int nb){
 double ans=0.;
 for(int e=0;e<nb;e++) Q->cnt[e].c=0;
 for(int e=0;e<N;e++) Q->cnt[b[e]-1].c++;
 for(int e=0;e<nb;e++){
  double pa=((double)Q->cnt[e].c)/((double)N);
  ans+=pa*pa;
 }
 return(ans);
}

