/* One mobile monomer on an open Lx x Ly board with optional removed (fixed) sites,
   interacting dimers with weight exp(-J N_parallel).
   Moves: (i) monomer slide m -> k through the dimer (j,k) covering a neighbour j of m,
          (ii) plaquette flip of two parallel dimers.  Slides alone are not ergodic
          (Bouttier, Bowick, Guitter, Jeng 2007); flips make the chain ergodic.
   usage: mono2 Lx Ly J nmoves seed nblocks out.bin [removed_site ...]
   output: nblocks x (Lx*Ly) doubles, monomer histogram per block (after 10% equilibration). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static int Lx,Ly,*mt,*rem; static uint64_t s0,s1;
static inline uint64_t rng(void){uint64_t a=s0,b=s1;s0=b;a^=a<<23;s1=a^b^(a>>17)^(b>>26);return s1+b;}
static inline double urand(void){return (rng()>>11)*(1.0/9007199254740992.0);}
static int nb(int i,int*o){int x=i%Lx,y=i/Lx,n=0;
  if(x+1<Lx&&!rem[i+1])o[n++]=i+1; if(x>0&&!rem[i-1])o[n++]=i-1;
  if(y+1<Ly&&!rem[i+Lx])o[n++]=i+Lx; if(y>0&&!rem[i-Lx])o[n++]=i-Lx; return n;}
static int plaq(int x,int y){int a=y*Lx+x,b=a+1,c=a+Lx,d=c+1;
  if(rem[a]||rem[b]||rem[c]||rem[d])return 0;
  return ((mt[a]==b&&mt[c]==d)||(mt[a]==c&&mt[b]==d))?1:0;}
static int local(int*s,int ns){int seen[128],n=0,tot=0;
  for(int u=0;u<ns;u++){int x=s[u]%Lx,y=s[u]/Lx;
    for(int dx=-1;dx<=0;dx++)for(int dy=-1;dy<=0;dy++){int px=x+dx,py=y+dy;
      if(px<0||py<0||px>Lx-2||py>Ly-2)continue;int id=py*Lx+px,dup=0;
      for(int v=0;v<n;v++)if(seen[v]==id)dup=1; if(!dup){seen[n++]=id;tot+=plaq(px,py);}}}
  return tot;}
/* greedy initial matching by augmenting paths (simple DFS), leaving one site free */
static int *vis;
static int aug(int u,int stamp){int o[4],n=nb(u,o);
  for(int t=0;t<n;t++){int v=o[t];if(vis[v]==stamp)continue;vis[v]=stamp;
    if(mt[v]<0||aug(mt[v],stamp)){mt[v]=u;mt[u]=v;return 1;}}
  return 0;}
int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage\n");return 1;}
  Lx=atoi(argv[1]);Ly=atoi(argv[2]);double J=atof(argv[3]);double nm=atof(argv[4]);
  int seed=atoi(argv[5]),nbk=atoi(argv[6]);char*out=argv[7];
  int N=Lx*Ly;mt=malloc(N*sizeof(int));rem=calloc(N,sizeof(int));vis=calloc(N,sizeof(int));
  for(int a=8;a<argc;a++)rem[atoi(argv[a])]=1;
  s0=0x9E3779B97F4A7C15ULL*(seed+1);s1=0xD1B54A32D192ED03ULL^s0;for(int i=0;i<20;i++)rng();
  for(int i=0;i<N;i++)mt[i]=-1;
  int stamp=0;
  for(int i=0;i<N;i++){int x=i%Lx,y=i/Lx;if(rem[i]||((x+y)&1))continue;if(mt[i]<0){stamp++;aug(i,stamp);}}
  int m=-1,nfree=0;for(int i=0;i<N;i++)if(!rem[i]&&mt[i]<0){m=i;nfree++;}
  if(nfree!=1){fprintf(stderr,"need exactly one free site, got %d\n",nfree);return 1;}
  double*hist=calloc((size_t)nbk*N,sizeof(double));
  long long neq=(long long)(0.1*nm),nmv=(long long)nm,accs=0,accf=0,ns=0,nf=0;
  long long per=nmv/nbk;
  for(long long t=0;t<neq+nmv;t++){
    if(urand()<0.5){ /* slide */
      ns++;int o[4],zm=nb(m,o);if(zm==0)goto rec;int j=o[(int)(urand()*zm)],kk=mt[j],o2[4],zk=nb(kk,o2);
      int s[3]={m,j,kk};int before=local(s,3);
      mt[m]=j;mt[j]=m;mt[kk]=-1;int after=local(s,3);
      double a=((double)zm/zk)*exp(-J*(after-before));
      if(a>=1||urand()<a){m=kk;accs++;}else{mt[kk]=j;mt[j]=kk;mt[m]=-1;}
    }else{ /* plaquette flip */
      nf++;int px=(int)(urand()*(Lx-1)),py=(int)(urand()*(Ly-1));
      int a=py*Lx+px,b=a+1,c=a+Lx,d=c+1;
      if(rem[a]||rem[b]||rem[c]||rem[d])goto rec;
      int hor=(mt[a]==b&&mt[c]==d),ver=(mt[a]==c&&mt[b]==d);
      if(!hor&&!ver)goto rec;
      int s[4]={a,b,c,d};int before=local(s,4);
      if(hor){mt[a]=c;mt[c]=a;mt[b]=d;mt[d]=b;}else{mt[a]=b;mt[b]=a;mt[c]=d;mt[d]=c;}
      int after=local(s,4);double w=exp(-J*(after-before));
      if(w>=1||urand()<w)accf++;
      else{if(hor){mt[a]=b;mt[b]=a;mt[c]=d;mt[d]=c;}else{mt[a]=c;mt[c]=a;mt[b]=d;mt[d]=b;}}
    }
    rec:
    if(t>=neq){long long bk=(t-neq)/per;if(bk>=nbk)bk=nbk-1;hist[bk*N+m]+=1;}
  }
  FILE*f=fopen(out,"wb");fwrite(hist,sizeof(double),(size_t)nbk*N,f);fclose(f);
  fprintf(stderr,"Lx=%d Ly=%d J=%g slide_acc=%.3f flip_acc=%.3f\n",Lx,Ly,J,(double)accs/ns,(double)accf/nf);
  return 0;}
