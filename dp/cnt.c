#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 I;
// usage: cnt R W r1 c1 r2 c2 ...   (row-major, mask over W columns)
int main(int argc,char**argv){int R=atoi(argv[1]),Wd=atoi(argv[2]),C=R*Wd,M=1<<Wd;
 int*blk=calloc(C,sizeof(int)); for(int i=3;i+1<argc;i+=2) blk[atoi(argv[i])*Wd+atoi(argv[i+1])]=1;
 I*F=calloc(M,sizeof(I)),*G=calloc(M,sizeof(I)); F[0]=1;
 for(int p=0;p<C;p++){int r=p/Wd,c=p%Wd; memset(G,0,M*sizeof(I));
  for(int m=0;m<M;m++){I v=F[m];if(!v)continue;
   if(blk[p]){if(!(m&1))G[m>>1]+=v;continue;}
   if(m&1){G[m>>1]+=v;continue;}
   if(r+1<R&&!blk[p+Wd])G[(m>>1)|(1<<(Wd-1))]+=v;
   if(c+1<Wd&&!blk[p+1]&&!(m&2))G[(m>>1)|1]+=v;}
  I*t=F;F=G;G=t;}
 printf("%.15Le\n",(long double)F[0]);}
