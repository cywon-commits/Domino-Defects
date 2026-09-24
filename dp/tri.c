#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// N x N triangular lattice (square grid + diagonal (r,c)-(r+1,c+1)), one monomer at v; count dimer covers.
// row-major; mask bit k (k=0..N) = cell p+k already covered.
typedef long double D;
int N; 
D count(int vr,int vc){int C=N*N,W=N+1,M=1<<W; D*F=calloc(M,sizeof(D)),*G=calloc(M,sizeof(D)); F[0]=1;
 for(int p=0;p<C;p++){int r=p/N,c=p%N; memset(G,0,M*sizeof(D));
  for(int m=0;m<M;m++){D v=F[m]; if(v==0) continue;
   int blk=(r==vr&&c==vc);
   if(m&1){ if(!blk) G[m>>1]+=v; continue;}
   if(blk){G[m>>1]+=v;continue;}
   int nm=m>>1;
   if(c+1<N && !(nm&1) && !(p+1==vr*N+vc)) G[nm|1]+=v;                 // right
   if(r+1<N && !(nm&(1<<(N-1))) && !(p+N==vr*N+vc)) G[nm|(1<<(N-1))]+=v; // up
   if(r+1<N && c+1<N && !(nm&(1<<N)) && !(p+N+1==vr*N+vc)) G[nm|(1<<N)]+=v; // diagonal
  }
  D*t=F;F=G;G=t;}
 D z=F[0]; free(F);free(G); return z;}
int main(int argc,char**argv){N=atoi(argv[1]); for(int r=0;r<N;r++)for(int c=0;c<N;c++) printf("%d %d %.12Le\n",r,c,count(r,c));}
