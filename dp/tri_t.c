#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Weighted version of tri.c: N x N square grid + diagonals (r,c)-(r+1,c+1) with weight t.
// Monomer at (vr,vc) (vr<0: no monomer). Prints Z. Reference for Pfaffian checks.
typedef long double D; int N; D T;
D count(int vr,int vc){int C=N*N,W=N+1,M=1<<W; D*F=calloc(M,sizeof(D)),*G=calloc(M,sizeof(D)); F[0]=1;
 int vp = vr<0 ? -1 : vr*N+vc;
 for(int p=0;p<C;p++){int r=p/N,c=p%N; memset(G,0,M*sizeof(D));
  for(int m=0;m<M;m++){D v=F[m]; if(v==0) continue;
   int blk=(p==vp);
   if(m&1){ if(!blk) G[m>>1]+=v; continue;}
   if(blk){G[m>>1]+=v;continue;}
   int nm=m>>1;
   if(c+1<N && !(nm&1) && p+1!=vp) G[nm|1]+=v;
   if(r+1<N && !(nm&(1<<(N-1))) && p+N!=vp) G[nm|(1<<(N-1))]+=v;
   if(T!=0 && r+1<N && c+1<N && !(nm&(1<<N)) && p+N+1!=vp) G[nm|(1<<N)]+=v*T;
  }
  D*t=F;F=G;G=t;}
 D z=F[0]; free(F);free(G); return z;}
int main(int argc,char**argv){N=atoi(argv[1]); T=strtold(argv[2],0);
 if(argc>4){printf("%.18Le\n",count(atoi(argv[3]),atoi(argv[4]))); return 0;}
 for(int r=0;r<N;r++)for(int c=0;c<N;c++) printf("%d %d %.18Le\n",r,c,count(r,c));}
