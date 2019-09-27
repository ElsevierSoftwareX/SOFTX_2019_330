struct ht3{
 struct ht3e *contents; //size N
 struct ht3e **map;
 int *counts;
};

struct ht3e{
 int x;
 int y;
 int z;
 int cXYZ;
 int *cXZ;
 int *cYZ;
 int *cZ;
 struct ht3e *nxt; //Next element of the same hash
};

struct ht3* R_allocHt3(int N){
 struct ht3 *ans=(struct ht3*)R_alloc(sizeof(struct ht3),1);
 ans->map=(struct ht3e**)R_alloc(sizeof(struct ht3e*),1024);
 ans->contents=(struct ht3e*)R_alloc(sizeof(struct ht3e),N);
 ans->counts=(int*)R_alloc(sizeof(int),N);
 return(ans);
}

int static inline disagree(struct ht3e *E,int x,int y,int z){
 return(!(
  ((x==0) || E->x==x) &&
  ((y==0) || E->y==y) &&
  ((z==0) || E->z==z)
 ));
}
int static inline hash(int x,int y,int z,int n){
 //return((x+y+z) & 0x3ff);
 //return(x & 0x3ff);
 //return((x+(y*4)+(z*8))%n);
 uint64_t pp=6364136223846793005;
 //uint64_t mix=((x*pp)+y)*pp+z;
 uint64_t mix=((x*pp)+x)*pp+x;
 uint32_t s=((mix>>18)^mix)>>27;
 uint32_t r=mix>>59;
 return ((s<<((-r)&31))|(s>>r))>>22;
}

//x &y start from 1
int fillHt3(struct ht3 *ht,int n,int nx,int *x,int ny,int *y,int nz,int *z){
 for(int e=0;e<1024;e++) ht->map[e]=NULL;
 int nE=0,nC=0;
 int *counts=ht->counts;

 for(int e=0;e<n;e++){
  int h=hash(x[e],y[e],z[e],n);
  struct ht3e **E;
  for(E=ht->map+h;
   (*E) &&
   x[e]!=(*E)->x &&
   y[e]!=(*E)->y &&
   z[e]!=(*E)->z;
   E=&((*E)->nxt)
  );

  if(!*E){
   //E not found, adding!
   *E=ht->contents+nE; nE++;
   (*E)->cXYZ=0;
   (*E)->nxt=NULL;
   (*E)->x=x[e];
   (*E)->y=y[e];
   (*E)->z=z[e];

   //Link XZ
   struct ht3e **EE;
   for(EE=ht->map+hash(x[e],0,z[e],n);(*EE)&&disagree(*EE,x[e],0,z[e]);EE=&((*EE)->nxt)); 
   if(!*EE || *EE==*E){
    *EE=*E;
    (*E)->cXZ=counts+nC; nC++;
    (*E)->cXZ[0]=0;
   }else (*E)->cXZ=(*EE)->cXZ;
 
   //Link YZ
   for(EE=ht->map+hash(0,y[e],z[e],n);(*EE)&&disagree(*EE,0,y[e],z[e]);EE=&((*EE)->nxt)); 
   if(!*EE || *EE==*E){
    *EE=*E;
    (*E)->cYZ=counts+nC; nC++;
    (*E)->cYZ[0]=0;
   }else (*E)->cYZ=(*EE)->cYZ;

   //Link Z
   for(EE=ht->map+hash(0,0,z[e],n);(*EE)&&disagree(*EE,0,0,z[e]);EE=&((*EE)->nxt)); 
   if(!*EE || *EE==*E){
    *EE=*E;
    (*E)->cZ=counts+nC; nC++;
    (*E)->cZ[0]=0;
   }else (*E)->cZ=(*EE)->cZ;
  }

  //E located, counting!
  (*E)->cXYZ++;
  (*E)->cXZ[0]++;
  (*E)->cYZ[0]++;
  (*E)->cZ[0]++;
 }
 return(nE);
}

double cmiHt3(struct ht3 *ht,int ne,int n){
 double I=0.,N=(double)n;
 for(int e=0;e<ne;e++){
  struct ht3e *E=ht->contents+e;
  double pXYZ=((double)E->cXYZ)/N;
  double Q=((double)E->cXYZ)*((double)E->cZ[0])/((double)E->cXZ[0])/((double)E->cYZ[0]);
  I+=pXYZ*log(Q);
 }
 return(I);
}
