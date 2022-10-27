#include <string.h>
#include <time.h>
#include <math.h>
#include <mygraph.h>
// Lattice Boltzmann for the diffusion equation in 1d

#define xdim 300
#define ydim 300
// velocities -1,0,1
int Xdim=xdim,Ydim=ydim;
typedef  double farr[3][3][xdim][ydim];  
farr *f,*ff,f1,f2;
double w[3][3]={{1./36.,1./9.,1./36.},{1./9.,4./9.,1./9.},{1./36.,1./9.,1./36.}},tau=1,rho[xdim][ydim],rhoth[xdim][ydim];
int t;

void getAnalytical(double t){
  double D=2./3*(tau-0.5);
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
    rhoth[x][y]=1+100/sqrt(2*3.14159*D*t)*exp(-((x-xdim/2)*(x-xdim/2))/(2*D*t))
      *1/sqrt(2*3.14159*D*t)*exp(-((y-ydim/2)*(y-ydim/2))/(2*D*t));
  }
}

void iterate(farr f,farr ff){
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      rho[x][y]=0;
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++){
	  rho[x][y]+=f[i][j][x][y];
	}
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++){
	  if ((pow(0.33*(x-xdim/2),2)/1600+pow(0.33*(y-ydim/2),2)/800*exp(0.01*0.33*(x-xdim/2)))<1) 
	    ff[i][j][(x+i-1+xdim)%xdim][(y+j-1+ydim)%ydim]=f[i][j][x][y]+
	      1/tau*(rho[x][y]*w[i][j]-f[i][j][x][y]);
	  else
	    ff[i][j][(x+i-1+xdim)%xdim][(y+j-1+ydim)%ydim]=100*w[i][j];
	}
    }
  t+=1;
}

void init(){
  t=0;
  f=&f1;
  ff=&f2;
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      rho[x][y]=100;
      if ((pow(0.33*(x-xdim/2),2)/1600+pow(0.33*(y-ydim/2),2)/800*exp(0.01*0.33*(x-xdim/2)))<1) rho[x][y]=20;
    }
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	  (*f)[i][j][x][y]=rho[x][y]*w[i][j];
  }
}

int main(){
    struct timespec TimeSpec;

    int sstep=0,cont=0,done=0,Repeat=1;
  TimeSpec.tv_sec=0;
  TimeSpec.tv_nsec=10000000;
  init();
  DefineGraphNxN_R("rho",&rho[0][0],&Xdim,&Ydim,NULL);
  DefineGraphNxN_R("rhoth",&rhoth[0][0],&Xdim,&Ydim,NULL);

  
  StartMenu("LB Diffusion",1);
  DefineInt("t",&t);
  DefineGraph(contour2d_,"Graph");
  DefineFunction("Reinitialize",init);
  DefineInt("Repeat",&Repeat);
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
      for(int i=0; i<Repeat; i++){
	sstep=0;
	iterate(*f,*ff);
	farr *tmp; // switching the f's
	tmp=f;
	f=ff;
	ff=tmp;
      }
    }
    
  }
}
