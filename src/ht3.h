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
 return((x+(y<<2)+(z<<3))%n);
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

   //Link XZ
   struct ht3e **EE;
   for(EE=ht->map+hash(x[e],0,z[e],n);(*EE)&&disagree(*EE,x[e],0,z[e]);EE=&((*EE)->nxt));
   if(!*EE){
    EE=E;
    (*E)->cXZ=counts+nC; nC++;
    (*E)->cXZ[0]=0;
   }else (*E)->cXZ=(*EE)->cXZ;
 
   //Link YZ
   for(EE=ht->map+hash(x[e],0,z[e],n);(*EE)&&disagree(*EE,x[e],0,z[e]);EE=&((*EE)->nxt));
   if(!*EE){
    EE=E;
    (*E)->cYZ=counts+nC; nC++;
    (*E)->cYZ[0]=0;
   }else (*E)->cYZ=(*EE)->cYZ;

   //Link Z
   for(EE=ht->map+hash(x[e],0,z[e],n);(*EE)&&disagree(*EE,x[e],0,z[e]);EE=&((*EE)->nxt));
   if(!*EE){
    EE=E;
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

void printHt3(struct ht3 *ht,int ne){
 for(int e=0;e<ne;e++){ 
  printf("State %d ",e);
  struct ht3e *E=ht->contents+e;
  printf("(%d %d %d) ",E->x,E->y,E->z);
  printf("cXYZ=%d ",E->cXYZ);
  printf("cXZ %p=%d ",E->cXZ,E->cXZ[0]);
  printf("cYZ %p=%d ",E->cYZ,E->cYZ[0]);
  printf("cZ %p=%d ",E->cZ,E->cZ[0]);
  printf("\n");
 }
}
