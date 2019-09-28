#define MAX(a,b) ((a)<(b))?(b):(a)
#define MIN(a,b) ((a)>(b))?(b):(a)

struct ht3{
 struct ht3e *contents; //size N
 struct ht3e **map;
 int *counts;
 int mBits;
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
 ans->contents=(struct ht3e*)R_alloc(sizeof(struct ht3e),N);
 ans->counts=(int*)R_alloc(sizeof(int),N*4);
 int mB=1;
 for(;(1<<mB)<=N;mB++);
 ans->mBits=MIN(14,mB+2);
 ans->map=(struct ht3e**)R_alloc(sizeof(struct ht3e*),1<<(ans->mBits));
 return(ans);
}

int hash(int x,int y,int z,int mB){
  return((x+y+z)%(1<<mB));
 uint64_t pp=6364136223846793005;
 uint64_t mix=((x*pp)+y)*pp+z;
 uint32_t s=((mix>>18)^mix)>>27;
 uint32_t r=mix>>59;
 return ((s<<((-r)&31))|(s>>r))%(1<<mB);//>>(32-mB);
}

//x &y start from 1
int fillHt3(struct ht3 *ht,int n,int nx,int *x,int ny,int *y,int nz,int *z){
 //Calculate maximal number of states
 uint64_t mStates,acc=((uint64_t)nx+1)*((uint64_t)ny+1);
 if(acc>0xffffffff){
  //Too many
  mStates=((uint64_t)1)<<63;
 }else{
  mStates=acc*((uint64_t)nz+1);
 }
 int mB=ht->mBits;
 if(mStates<(1<<mB)){
  printf("Shortcut possible :(\n");
 }
 for(int e=0;e<(1<<mB);e++) ht->map[e]=NULL;
 int nE=0,nC=0;
 int *counts=ht->counts;

 for(int e=0;e<n;e++){
  int X=x[e],Y=y[e],Z=z[e];
  struct ht3e **E;
  for(E=ht->map+hash(X,Y,Z,mB);(*E) && (X!=(*E)->x || Y!=(*E)->y || Z!=(*E)->z);E=&((*E)->nxt));

  if(!*E){
   //E not found, adding!
   *E=ht->contents+nE; nE++;
   (*E)->cXYZ=0;
   (*E)->nxt=NULL;
   (*E)->x=X; (*E)->y=Y; (*E)->z=Z;

   //Link XZ
   struct ht3e **EE;
   for(EE=ht->map+hash(X,0,Z,mB);(*EE) && (X!=(*EE)->x || Z!=(*EE)->z);EE=&((*EE)->nxt));

   if(!*EE || *EE==*E){
    *EE=*E;
    (*E)->cXZ=counts+nC; nC++;
    (*E)->cXZ[0]=0;
   }else (*E)->cXZ=(*EE)->cXZ;
 
   //Link YZ
   for(EE=ht->map+hash(0,Y,Z,mB);(*EE) && (Y!=(*EE)->y || Z!=(*EE)->z);EE=&((*EE)->nxt));
   if(!*EE || *EE==*E){
    *EE=*E;
    (*E)->cYZ=counts+nC; nC++;
    (*E)->cYZ[0]=0;
   }else (*E)->cYZ=(*EE)->cYZ;

   //Link Z
   for(EE=ht->map+hash(0,0,Z,mB);(*EE) && Z!=(*EE)->z;EE=&((*EE)->nxt));
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
