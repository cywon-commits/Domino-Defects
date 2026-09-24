#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 I;
int main(int argc,char**argv){int N=atoi(argv[1]),C=N*N,M=1<<N;
 I*f=calloc((size_t)(C+1)*M,sizeof(I)),*g=calloc((size_t)(C+1)*M,sizeof(I));int*blk=calloc(C,sizeof(int));blk[0]=1;
 f[0]=1;
 for(int p=0;p<C;p++){int r=p/N,c=p%N;I*F=f+(size_t)p*M,*G=f+(size_t)(p+1)*M;
  for(int m=0;m<M;m++){I v=F[m];if(!v)continue;
   if(blk[p]){if(!(m&1))G[m>>1]+=v;continue;}
   if(m&1){G[m>>1]+=v;continue;}
   if(r+1<N&&!blk[p+N])G[(m>>1)|(1<<(N-1))]+=v;
   if(c+1<N&&!blk[p+1]&&!(m&2))G[(m>>1)|1]+=v;}}
 g[(size_t)C*M]=1;
 for(int p=C-1;p>=0;p--){int r=p/N,c=p%N;I*G=g+(size_t)p*M,*H=g+(size_t)(p+1)*M;
  for(int m=0;m<M;m++){I v=0;
   if(blk[p]){if(!(m&1))v=H[m>>1];G[m]=v;continue;}
   if(m&1){G[m]=H[m>>1];continue;}
   if(r+1<N&&!blk[p+N])v+=H[(m>>1)|(1<<(N-1))];
   if(c+1<N&&!blk[p+1]&&!(m&2))v+=H[(m>>1)|1];G[m]=v;}}
 for(int b=1;b<C;b++){if(((b/N+b%N)&1)==0)continue;I s=0;I*F=f+(size_t)b*M,*H=g+(size_t)(b+1)*M;
  for(int m=0;m<M;m+=2)s+=F[m]*H[m>>1];
  unsigned long long hi=(unsigned long long)(s>>64), lo=(unsigned long long)s; printf("%d %.12Le %llu %llu\n",b,(long double)s,hi,lo);}
}
