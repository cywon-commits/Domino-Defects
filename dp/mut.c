#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 I;
int N; int *blk;
int main(int argc,char**argv){
  N=atoi(argv[1]); int C=N*N, M=1<<N;
  I *f=malloc(sizeof(I)*(C+1)*M), *g=malloc(sizeof(I)*(C+1)*M);
  blk=calloc(C,sizeof(int));
  long double *W=calloc((size_t)C*C,sizeof(long double));
  for(int a=0;a<C;a++){ int ar=a/N,ac=a%N; if(((ar+ac)&1)==0) continue;
    memset(blk,0,C*sizeof(int)); blk[0]=1; blk[C-1]=1; blk[a]=1;
    memset(f,0,sizeof(I)*(C+1)*M); memset(g,0,sizeof(I)*(C+1)*M);
    f[0]=1;
    for(int p=0;p<C;p++){int r=p/N,c=p%N; I*F=f+(size_t)p*M,*G=f+(size_t)(p+1)*M;
      for(int m=0;m<M;m++){I v=F[m]; if(!v) continue;
        if(blk[p]){ if(!(m&1)) G[m>>1]+=v; continue;}
        if(m&1){G[m>>1]+=v;continue;}
        if(r+1<N && !blk[p+N]) G[(m>>1)|(1<<(N-1))]+=v;
        if(c+1<N && !blk[p+1] && !(m&2)) G[(m>>1)|1]+=v; }}
    g[(size_t)C*M+0]=1;
    for(int p=C-1;p>=0;p--){int r=p/N,c=p%N; I*G=g+(size_t)p*M,*H=g+(size_t)(p+1)*M;
      for(int m=0;m<M;m++){I v=0;
        if(blk[p]){ if(!(m&1)) v=H[m>>1]; G[m]=v; continue;}
        if(m&1){G[m]=H[m>>1];continue;}
        if(r+1<N && !blk[p+N]) v+=H[(m>>1)|(1<<(N-1))];
        if(c+1<N && !blk[p+1] && !(m&2)) v+=H[(m>>1)|1];
        G[m]=v;}}
    for(int b=a+1;b<C;b++){int br=b/N,bc=b%N; if(((br+bc)&1)==0) continue;
      I s=0; I*F=f+(size_t)b*M,*H=g+(size_t)(b+1)*M;
      for(int m=0;m<M;m+=2) s+=F[m]*H[m>>1];
      W[(size_t)a*C+b]=(long double)s; }
  }
  for(int a=0;a<C;a++)for(int b=a+1;b<C;b++) if(W[(size_t)a*C+b]>0) printf("%d %d %.10Le\n",a,b,W[(size_t)a*C+b]);
}
