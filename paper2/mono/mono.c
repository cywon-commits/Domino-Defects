/* One forced monomer on an open odd L x L board, interacting dimers exp(-J N_parallel).
   usage: mono L J nmoves seed out.bin */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static int L,*mt; static uint64_t s0,s1;
static inline uint64_t rng(void){uint64_t a=s0,b=s1;s0=b;a^=a<<23;s1=a^b^(a>>17)^(b>>26);return s1+b;}
static inline double urand(void){return (rng()>>11)*(1.0/9007199254740992.0);}
static int nb(int i,int*o){int x=i%L,y=i/L,n=0;if(x+1<L)o[n++]=i+1;if(x>0)o[n++]=i-1;if(y+1<L)o[n++]=i+L;if(y>0)o[n++]=i-L;return n;}
static int plaq(int x,int y){int a=y*L+x,b=a+1,c=a+L,d=c+1;return ((mt[a]==b&&mt[c]==d)||(mt[a]==c&&mt[b]==d))?1:0;}
static int local(int*s,int ns){ /* plaquettes touching any of the sites, deduplicated */
  int seen[64],n=0,tot=0;
  for(int u=0;u<ns;u++){int x=s[u]%L,y=s[u]/L;
    for(int dx=-1;dx<=0;dx++)for(int dy=-1;dy<=0;dy++){int px=x+dx,py=y+dy;if(px<0||py<0||px>L-2||py>L-2)continue;
      int id=py*L+px,dup=0;for(int v=0;v<n;v++)if(seen[v]==id)dup=1;if(!dup){seen[n++]=id;tot+=plaq(px,py);}}}
  return tot;}
int main(int argc,char**argv){
  L=atoi(argv[1]);double J=atof(argv[2]);double nm=atof(argv[3]);
  s0=0x9E3779B97F4A7C15ULL*(atoi(argv[4])+1);s1=0xD1B54A32D192ED03ULL^s0;for(int i=0;i<20;i++)rng();
  int N=L*L;mt=malloc(N*sizeof(int));double*hist=calloc(N,sizeof(double));
  int path[N],k=0;for(int y=0;y<L;y++){if(y%2==0)for(int x=0;x<L;x++)path[k++]=y*L+x;else for(int x=L-1;x>=0;x--)path[k++]=y*L+x;}
  for(int i=0;i+1<N;i+=2){mt[path[i]]=path[i+1];mt[path[i+1]]=path[i];}
  int m=path[N-1];mt[m]=-1;
  long long neq=(long long)(0.1*nm),nmv=(long long)nm,acc=0;
  for(long long t=0;t<neq+nmv;t++){
    int o[4],zm=nb(m,o),j=o[(int)(urand()*zm)],kk=mt[j],o2[4],zk=nb(kk,o2);
    int s[3]={m,j,kk};int before=local(s,3);
    mt[m]=j;mt[j]=m;mt[kk]=-1;int after=local(s,3);
    double a=((double)zm/zk)*exp(-J*(after-before));
    if(a>=1||urand()<a){m=kk;acc++;}else{mt[kk]=j;mt[j]=kk;mt[m]=-1;}
    if(t>=neq)hist[m]+=1;
  }
  FILE*f=fopen(argv[5],"wb");fwrite(hist,sizeof(double),N,f);fclose(f);
  fprintf(stderr,"L=%d J=%g acc=%.3f\n",L,J,(double)acc/(neq+nmv));return 0;}
