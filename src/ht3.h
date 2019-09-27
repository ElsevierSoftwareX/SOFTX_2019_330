struct ht3{
 struct ht3e *contents; //size N
 struct ht3e **map;
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

int static inline disagree(struct ht3e *E,int x,int y,int z){
 return(!(
  ((x==0) || E->x==x) &&
  ((y==0) || E->y==y) &&
  ((z==0) || E->z==z)
 ));
}
int static inline hash(int x,int y,int z,int n){
 return((x+(y*4)+(z*8))%n);
}

//x &y start from 1
int fillHt3(struct ht3 *ht,int n,int nx,int *x,int ny,int *y,int nz,int *z,int *counts){
 for(int e=0;e<n;e++) ht->map[e]=NULL;
 int nE=0,nC=0;

 for(int e=0;e<n;e++){
  int h=hash(x[e],y[e],z[e],n);
  struct ht3e **E;
  for(E=ht->map+h;(*E)&&disagree(*E,x[e],y[e],z[e]);E=&((*E)->nxt));

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
