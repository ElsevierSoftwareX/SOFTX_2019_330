uint32_t rng(uint64_t *state,uint32_t stream){
 uint32_t rot=state[0]>>59;
 uint32_t s=((state[0]>>18)^state[0])>>27;
 state[0]=state[0]*6364136223846793005+(2*stream+1);
 return (s<<((-rot)&31))|(s>>rot);
}
