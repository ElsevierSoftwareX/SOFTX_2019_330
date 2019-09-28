#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//Shims
char* R_alloc(size_t a,size_t b){
 return(malloc(a*b));
}
double unif_rand(){
 return(.7);
}
int min(int a,int b){
 return(a<b?a:b);
}
void GetRNGstate(){}
void PutRNGstate(){}
#define NA_INTEGER 0xffffffff
void error(const char *e){
 abort();
}

#include "../src/ht3.h"
#include "../src/ht.h"
#include "../src/pcg.h"

int main(int argc,char **argv){
 printf("mb\tn\tnx\tny\tnz\tns\tcmi\tnst\toccp\tpoccp\tmD\taD\ttime\ttimeOld\n");
 uint64_t r=21; rng(&r,0);
 struct timespec tpa,tpb;
 double timeOld=0.;
 double timeNew=0.;

 for(int j=0;j<1;j++){

  
  int n,nx,ny,nz;
  if(1){
   if(rng(&r,0)%2==0) n=50+rng(&r,0)%100; else n=100+rng(&r,0)%300;
   nx=min(1+rng(&r,0)%50,n);
   ny=min(1+rng(&r,0)%50,n);
   nz=min(1+rng(&r,0)%50,n);
  }else{
   //Madelon
   n=2000;
   nx=10;
   ny=10;
   nz=2;
  }
  int *x=malloc(sizeof(int)*n),
   *y=malloc(sizeof(int)*n),
   *z=malloc(sizeof(int)*n);
  struct ht3 *H=R_allocHt3(n);

  struct ht *h=R_allocHt(n);
  int *xz=malloc(sizeof(int)*n);
  int *cXZ=malloc(sizeof(int)*n);
  int *cY=malloc(sizeof(int)*n);
  
  for(int e=0;e<n;e++){
   x[e]=rng(&r,0)%nx+1;
   y[e]=rng(&r,0)%ny+1;
   z[e]=rng(&r,0)%nz+1;
   int nB=5;
   printf(">>> %d %d %d => %d %d %d %d\n",x[e],y[e],z[e],hash(x[e],y[e],z[e],nB),hash(x[e],0,z[e],nB),hash(0,y[e],z[e],nB),hash(z[e],1,z[e],nB));
  }
 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tpb);
  double timeHt3=0.; int ni=0;
  int ns; double cmi;
  while(timeHt3<10000){
   ni++;
   ns=fillHt3(H,n,nx,x,ny,y,nz,z);
   cmi=cmiHt3(H,ns,n);
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tpa);
   timeHt3=((tpa.tv_sec-tpb.tv_sec)*1e9+(tpa.tv_nsec-tpb.tv_nsec));
  }
  timeHt3/=((double)ni)*((double)n);
  int numOccupied=0;
  int maxDepth=0;
  double aveDepth=0;
  for(int e=0;e<(1<<H->mBits);e++){
   numOccupied+=!!H->map[e];
   int depth=0;
   for(struct ht3e *a=H->map[e];a;a=a->nxt) depth++;
   aveDepth+=depth;
   maxDepth=MAX(maxDepth,depth);
  }
  double occupation=((double)numOccupied)/((double)(1<<H->mBits));
  aveDepth/=numOccupied;

  //Alternative ht
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tpb);
  double timeHt=0.; ni=0;
  while(timeHt<10000){
   ni++;
   int nxz=fillHt(h,n,nz,z,nx,x,xz,NULL,NULL,1);
   fillHt(h,n,ny,y,nxz,xz,NULL,cY,cXZ,0); //filling cY is a bit cheating
   miHt(h,cY,cXZ);
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tpa);
   timeHt=((tpa.tv_sec-tpb.tv_sec)*1e9+(tpa.tv_nsec-tpb.tv_nsec));
  }
  timeHt/=((double)ni)*((double)n);

  timeNew+=timeHt3;
  timeOld+=timeHt;


  printf("%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%0.3f\t%d\t%0.3f\t%0.6f\t%0.6f\n",H->mBits,n,nx,ny,nz,ns,cmi,ns,numOccupied,occupation,maxDepth,aveDepth,timeHt3,timeHt);

  free(x);
  free(y);
  free(z);
  free(H->counts);
  free(H->contents);
  free(H->map);
  free(H);

  free(xz);
  free(h->map);
  free(h->cnt);
  free(h);
 }
 fprintf(stderr,"Slowdown: %0.3f\n",timeNew/timeOld);
 
 return(0);
}
