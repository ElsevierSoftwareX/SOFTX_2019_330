#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

#include "../src/ht3.h"
#include "../src/pcg.h"

int main(int argc,char **argv){
 uint64_t r=21; rng(&r,0);

 int n=10+rng(&r,0)%30;
 int nx=min(1+rng(&r,0)%20,n);
 int ny=min(1+rng(&r,0)%20,n);
 int nz=min(1+rng(&r,0)%20,n);
 int *x=malloc(sizeof(int)*n),
  *y=malloc(sizeof(int)*n),
  *z=malloc(sizeof(int)*n);
 struct ht3 *H=R_allocHt3(n);
 
 for(int e=0;e<n;e++){
  x[e]=rng(&r,0)%nx+1;
  y[e]=rng(&r,0)%ny+1;
  z[e]=rng(&r,0)%nz+1;
 }

 int ns=fillHt3(H,n,nx,x,ny,y,nz,z);
 printf("CMI[%d]=%0.3f\n",ns,cmiHt3(H,ns,n));
 
 return(0);
}
