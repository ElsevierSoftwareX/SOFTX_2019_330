uint64_t seed_from_r(){
 GetRNGstate();
 uint64_t a=((uint64_t)(((double)(~((uint32_t)0)))*unif_rand()))<<32;
 a+=(uint32_t)(((double)(~((uint32_t)0)))*unif_rand());
 PutRNGstate();
 return(a);
}

uint32_t rng(uint64_t *state,uint32_t stream){
 state[0]=state[0]*6364136223846793005+(2*stream+1);
 uint32_t rot=state[0]>>59;
 uint32_t s=((state[0]>>18)^state[0])>>27;
 return (s<<((-rot)&31))|(s>>rot);
}
