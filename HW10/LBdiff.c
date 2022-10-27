#include <string.h>
#include <time.h>
#include <math.h>
#include <mygraph.h>
// Lattice Boltzmann for the diffusion equation in 1d

#define xdim 100
// velocities -1,0,1
int Xdim=xdim;
double f[3][xdim],w[3]={1./6.,2./3.,1./6.},tau=1,rho[xdim],rhoth[xdim],t;

void getAnalytical(double t){
  double D=2./3*(tau-0.5);
  for (int x=0; x<xdim; x++){
    rhoth[x]=1+100/sqrt(2*3.14159*D*t)*exp(-((x-xdim/2)*(x-xdim/2))/(2*D*t));
  }
}

void iterate(double f[3][xdim]){
  for (int x=0; x<xdim; x++){
    rho[x]=f[0][x]+f[1][x]+f[2][x];
    for (int v=0; v<3; v++) f[v][x]+=1/tau*(rho[x]*w[v]-f[v][x]);
  }
  double tmp=f[0][0];
  memmove(&f[0][0],&f[0][1],(xdim-1)*sizeof(double));
  f[0][xdim-1]=tmp; // for periodic boundary conditions.
  tmp=f[2][xdim-1];
  memmove(&f[2][1],&f[2][0],(xdim-1)*sizeof(double));
  f[2][0]=tmp; // for periodic boundary conditions.
  t+=1;
}

void init(){
  t=0;
  for (int x=0; x<xdim; x++){
    rho[x]=1;
    if (x==xdim/2) rho[x]=100;
  }
  for (int x=0; x<xdim; x++){
    for (int v=0; v<3; v++) f[v][x]=rho[x]*w[v];
  }
}

int main(){
    struct timespec TimeSpec;

  int sstep=0,cont=0,done=0;
  TimeSpec.tv_sec=0;
  TimeSpec.tv_nsec=10000000;
  init();
  DefineGraphN_R("rho",&rho[0],&Xdim,NULL);
  DefineGraphN_R("rhoth",&rhoth[0],&Xdim,NULL);

  
  StartMenu("LB Diffusion",1);
  DefineGraph(curve2d_,"Graph");
  DefineFunction("Reinitialize",init);
  DefineDouble("tau",&tau);
  DefineBool("sstep",&sstep);
  DefineBool("cont",&cont);
  DefineBool("done",&done);
  EndMenu();

  while (!done){
    getAnalytical(t);
    Events(1);
    DrawGraphs();
    nanosleep(&TimeSpec,NULL);
 
    if (sstep||cont){
      sstep=0;
      iterate(f);
    }
    
  }
}
