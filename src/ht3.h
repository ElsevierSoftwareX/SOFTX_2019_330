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
 struct ht3e *nxt; //Next element of the same hash
};

struct ht3* R_allocHt3(int N){
 struct ht3 *ans=(struct ht3*)R_alloc(sizeof(struct ht3),1);
 ans->contents=(struct ht3e*)R_alloc(sizeof(struct ht3e),N);
 ans->counts=(int*)R_alloc(sizeof(int),N*4);
 int mB=1;
 for(;(1<<mB)<=N;mB++);
 ans->mBits=MIN(17,mB+3);
 ans->map=(struct ht3e**)R_alloc(sizeof(struct ht3e*),1<<(ans->mBits));
 return(ans);
}

int hash(int x,int y,int z,int mB){
 uint64_t pp=6364136223846793005;
 uint64_t mix=((z*pp+y)*pp+x)*pp;
 return mix%(1<<mB);
}
#define PH(x,y,z) ((x)*(ny+1)*(nz)+(y)*(nz)+(z-1))

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
  //Shortcut, we ignore hash part
  for(int e=0;e<mStates;e++) ht->map[e]=NULL;
  int nE=0,nC=0;
  int *counts=ht->counts;

  for(int e=0;e<n;e++){
   int X=x[e],Y=y[e],Z=z[e];
   if(PH(X,Y,Z)>mStates) printf("EEEE!");
   struct ht3e **E=ht->map+PH(X,Y,Z);
   if(!*E){
    //E not found, adding!
    *E=ht->contents+nE; nE++;
    (*E)->cXYZ=0;
    (*E)->x=X; (*E)->y=Y; (*E)->z=Z;

    //Link XZ
    struct ht3e **EE;
    EE=ht->map+PH(X,0,Z);

    if(!*EE){
     *EE=*E;
     (*E)->cXZ=counts+nC; nC++;
     (*E)->cXZ[0]=0;
    }else (*E)->cXZ=(*EE)->cXZ;
  
    //Link YZ
    EE=ht->map+PH(0,Y,Z);
    if(!*EE){
     *EE=*E;
     (*E)->cYZ=counts+nC; nC++;
     (*E)->cYZ[0]=0;
    }else (*E)->cYZ=(*EE)->cYZ;
   }

   //E located, counting!
   (*E)->cXYZ++;
   (*E)->cXZ[0]++;
   (*E)->cYZ[0]++;
  }
  return(nE);
 }else{
  //Real HT
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
   }

   //E located, counting!
   (*E)->cXYZ++;
   (*E)->cXZ[0]++;
   (*E)->cYZ[0]++;
  }
  return(nE);
 }
}

double cmiHt3(struct ht3 *ht,int ne,int n,int *cZ){
 double I=0.,N=(double)n;
 for(int e=0;e<ne;e++){
  struct ht3e *E=ht->contents+e;
  double pXYZ=((double)E->cXYZ)/N;
  double Q=((double)E->cXYZ)*((double)cZ[E->z-1])/((double)E->cXZ[0])/((double)E->cYZ[0]);
  I+=pXYZ*log(Q);
 }
 return(I);
}
