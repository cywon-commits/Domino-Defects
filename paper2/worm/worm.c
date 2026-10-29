/* Interacting dimers (Alet et al.): weight exp(-J * N_parallel_plaquettes), J = beta*v.
   Cylinder: W periodic (x), H open (y). Two monomers: tail fixed at (xt,yt), head moves.
   Head move: pick neighbour j of head i uniformly (z_i), k = partner(j); rematch i-j, head -> k.
   Acceptance min(1, z_i/z_k * exp(-J dN)). Histogram of head position = Z(tail, head) up to a constant.
   usage: worm W H J xt yt nmoves seed out.bin */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static int W,H,*mt; static uint64_t s0,s1;
static inline uint64_t rng(void){uint64_t a=s0,b=s1;s0=b;a^=a<<23;s1=a^b^(a>>17)^(b>>26);return s1+b;}
static inline double urand(void){return (rng()>>11)*(1.0/9007199254740992.0);}
static inline int id(int x,int y){return y*W+((x%W)+W)%W;}
static int nb(int i,int*o){int x=i%W,y=i/W,n=0;o[n++]=id(x+1,y);o[n++]=id(x-1,y);if(y+1<H)o[n++]=id(x,y+1);if(y>0)o[n++]=id(x,y-1);return n;}
static int plaq(int x,int y){ /* plaquette with lower-left corner (x,y), 0<=y<=H-2 */
  int a=id(x,y),b=id(x+1,y),c=id(x,y+1),d=id(x+1,y+1);
  return ((mt[a]==b&&mt[c]==d)||(mt[a]==c&&mt[b]==d))?1:0;}
static int bondplaq(int a,int b,int*p){ /* plaquette ids (x + W*y) touching bond a-b */
  int n=0,xa=a%W,ya=a/W,xb=b%W,yb=b/W;
  if(ya==yb){int xl=(xb==(xa+1)%W)?xa:xb; if(ya<=H-2)p[n++]=ya*W+xl; if(ya>=1)p[n++]=(ya-1)*W+xl;}
  else{int yl=ya<yb?ya:yb; p[n++]=yl*W+xa; p[n++]=yl*W+((xa-1+W)%W);}
  return n;}
static int sumplaq(int*p,int n){int s=0;for(int u=0;u<n;u++){int dup=0;for(int v=0;v<u;v++)if(p[v]==p[u])dup=1;if(!dup)s+=plaq(p[u]%W,p[u]/W);}return s;}
int main(int argc,char**argv){
  W=atoi(argv[1]);H=atoi(argv[2]);double J=atof(argv[3]);int xt=atoi(argv[4]),yt=atoi(argv[5]);
  double nm=atof(argv[6]);s0=0x9E3779B97F4A7C15ULL*(atoi(argv[7])+1);s1=0xD1B54A32D192ED03ULL^s0;for(int i=0;i<20;i++)rng();
  int N=W*H; mt=malloc(N*sizeof(int)); double*hist=calloc(N,sizeof(double));
  for(int y=0;y<H;y++)for(int x=0;x<W;x+=2){mt[id(x,y)]=id(x+1,y);mt[id(x+1,y)]=id(x,y);}
  int tail=id(xt,yt),head=mt[tail]; mt[tail]=-1; mt[head]=-1;
  double bf[9]; for(int d=-4;d<=4;d++) bf[d+4]=exp(-J*d);
  long long acc=0,tot=0; long long neq=(long long)(nm*0.1), nmv=(long long)nm;
  for(long long t=0;t<neq+nmv;t++){
    int o[4],zi=nb(head,o),j=o[(int)(urand()*zi)];
    if(j!=tail){
      int k=mt[j],o2[4],zk=nb(k,o2),p[8],np=0;
      np+=bondplaq(j,k,p+np); np+=bondplaq(head,j,p+np);
      int before=sumplaq(p,np);
      mt[head]=j;mt[j]=head;mt[k]=-1;
      int after=sumplaq(p,np);
      double a=((double)zi/zk)*bf[after-before+4];
      if(a>=1.0||urand()<a){head=k;acc++;}
      else{mt[k]=j;mt[j]=k;mt[head]=-1;}
    }
    if(t>=neq){hist[head]+=1.0;tot++;}
  }
  FILE*f=fopen(argv[8],"wb");fwrite(hist,sizeof(double),N,f);fclose(f);
  fprintf(stderr,"W=%d H=%d J=%g tail=(%d,%d) acc=%.3f\n",W,H,J,xt,yt,(double)acc/(neq+nmv));
  return 0;}
