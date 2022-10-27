#include <string.h>
#include <time.h>
#include <math.h>
#include <mygraph.h>
// Lattice Boltzmann for the diffusion equation in 1d

#define xdim 100
#define ydim 100
#define Links 10000
// velocities -1,0,1
int Xdim=xdim,Ydim=ydim;
typedef  double farr[3][3][xdim][ydim];  
farr *f,*ff,f1,f2,*g,*gg,g1,g2;
//Phi is the density of gas
double w[3][3]={{1./36.,1./9.,1./36.},{1./9.,4./9.,1./9.},{1./36.,1./9.,1./36.}},tau=1,tauphiin=1,tauphiout=0.6,rho[xdim][ydim],phi[xdim][ydim],u[xdim][ydim][2],rhosym[xdim][ydim],usym[xdim][ydim][2],t,Acc=0,hole=5,heat=0,cool=0;
int rhosymreq=0,usymreq=0;
int v[3]={-1,0,1};

// Library for broken links. Make sure you only count them once!
int link[Links][3],links=0,field[xdim][ydim];
double fieldg[xdim][ydim];
int fieldgreq=0;

//This defines several "blocks" where I can do different things
void DefineField(){
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      field[x][y]=1;
      if (y-(ydim/10)<0) field[x][y]=2;
      if (y-(9*ydim/10)>0) field[x][y]=3;
      //if ((y-ydim/2)*(y-ydim/2)<200-80*cos(2*3.14159*x/xdim)) field[x][y]=1;
      // if ((x-xdim/2)*(x-xdim/2)/4+(y-ydim/2)*(y-ydim/2)<100) field[x][y]=2;
    }
}
/*
int boundary(int x1,int y1,int x2, int y2){
  int xc=xdim/2; int yc=ydim/2;
  int r2=(xdim/3)*(xdim/3);
  if (((x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-r2)*((x2-xc)*(x2-xc)+(y2-yc)*(y2-yc)-r2)<=0){
    if (((x1-xc)<(y1-yc+hole))&&((x1-yc)>(y1-yc-hole))
	&&((x2-xc)<(y2-yc+hole))&&((x2-yc)>(y2-yc-hole)))
      return 0;
    else
      return 1;
  }
  else
    return 0;
}
*/
int boundary(int x1,int y1, int x2, int y2){
  if (field[x1][y1]!=field[x2][y2]) return 1;
  else return 0;
}

void defineBoundary(){
  links=0;
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      if (boundary(x,y,(x+1)%xdim,y)){
	link[links][0]=0;
	link[links][1]=x;
	link[links][2]=y;
	links++;
      }
      if (boundary(x,y,(x+1)%xdim,(y+ydim-1)%ydim)){
	link[links][0]=1;
	link[links][1]=x;
	link[links][2]=y;
	links++;
      }
      if (boundary(x,y,x,(y+1)%ydim)){
	link[links][0]=2;
	link[links][1]=x;
	link[links][2]=y;
	links++;
      }
      if (boundary(x,y,(x+1)%xdim,(y+1)%ydim)){
	link[links][0]=3;
	link[links][1]=x;
	link[links][2]=y;
	links++;
      }
    }
  if (links>Links) printf("Warning increase Links!!!");
}

void bounceback(farr f){
  double tmpf;
  for (int l=0; l<links; l++){
    int x=link[l][1];
    int y=link[l][2];
    int xp=(x+1)%xdim;
    int yp=(y+1)%ydim;
    int ym=(y+ydim-1)%ydim;
    switch (link[l][0]){
    case 0:// 1 0
      tmpf=f[0][1][x][y];f[0][1][x][y]=f[2][1][xp][y]; f[2][1][xp][y]=tmpf; 
      break;
    case 1:// 1 -1
      tmpf=f[0][2][x][y];f[0][2][x][y]=f[2][0][xp][ym];f[2][0][xp][ym]=tmpf;
      break;
    case 2:// 0 1
      tmpf=f[1][0][x][y];f[1][0][x][y]=f[1][2][x][yp];f[1][2][x][yp]=tmpf;
      break;
    case 3:// 1 1
      tmpf=f[0][0][x][y];f[0][0][x][y]=f[2][2][xp][yp];f[2][2][xp][yp]=tmpf;
      break;
    }
  }
}

void gboundary(farr g){
  for(int x=0; x<xdim; x++){
    double tmp=g[0][2][x][0];
    g[0][2][x][0]=g[2][0][(x+1)%xdim][ydim-1];
    g[2][0][(x+1)%xdim][ydim-1]=tmp;

    tmp=g[1][2][x][0];
    g[1][2][x][0]=g[1][0][x][ydim-1];
    g[1][0][x][ydim-1]=tmp;

    tmp=g[2][2][x][0];
    g[2][2][x][0]=g[0][0][(x-1+xdim)%xdim][ydim-1];
    g[0][0][(x-1+xdim)%xdim][ydim-1]=tmp;
  }
}

double feq(int i,int j, double rho, double ux, double uy){
  return rho*w[i][j]*(1+3*(v[i]*ux+v[j]*uy)+9./2*(v[i]*ux+v[j]*uy)*(v[i]*ux+v[j]*uy)-3./2*(ux*ux+uy*uy));
}


void iterate(farr f,farr ff,farr g, farr gg){
  double tauphi;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int x=0; x<xdim; x++){
	g[i][j][x][0]+=feq(i,j,heat,0,0);
	g[i][j][x][ydim-1]=feq(i,j,cool,0,0);
      }
    }
  }
  
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      rho[x][y]=0;
      phi[x][y]=0;
      u[x][y][0]=0;
      u[x][y][1]=0;
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++){
	  rho[x][y]+=f[i][j][x][y];
	  phi[x][y]+=g[i][j][x][y];
	  u[x][y][0]+=f[i][j][x][y]*v[i];
	  u[x][y][1]+=f[i][j][x][y]*v[j];
	}
      u[x][y][0]/=rho[x][y];
      u[x][y][1]/=rho[x][y];
      if (field[x][y]!=1){
	u[x][y][0]=0;
	u[x][y][1]=0;
	tauphi=tauphiin;
      }
      else tauphi=tauphiout;
      
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++){
	  // Forcing term with constant force (like gravity)
	  f[i][j][x][y]+=feq(i,j,rho[x][y],u[x][y][0],u[x][y][1]+Acc*(phi[x][y]-0.75))
	    -feq(i,j,rho[x][y],u[x][y][0],u[x][y][1]);
	  //g[i][j][x][y]+=feq(i,j,phi[x][y],u[x][y][0],u[x][y][1]+Acc*(phi[x][y]-0.75))-feq(i,j,phi[x][y],u[x][y][0],u[x][y][1]);
	  ff[i][j][(x+i-1+xdim)%xdim][(y+j-1+ydim)%ydim]=f[i][j][x][y]+
			     1/tau*(feq(i,j,rho[x][y],u[x][y][0],u[x][y][1]+Acc*(phi[x][y]-0.75))-f[i][j][x][y]);
	  gg[i][j][(x+i-1+xdim)%xdim][(y+j-1+ydim)%ydim]=g[i][j][x][y]+
			     1/tauphi*(feq(i,j,phi[x][y],u[x][y][0],u[x][y][1]+Acc*(phi[x][y]-0.75))-g[i][j][x][y]);
	}
      //if(field[x][y]==2) phi[x][y]=2;
      //if(field[x][y]==3) phi[x][y]=0.5;
    }
  bounceback(ff);
  //bounceback(gg);
  gboundary(gg);
  t+=1;
}

void init(){
  t=0;
  f=&f1;
  ff=&f2;
  g=&g1;
  gg=&g2;
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      rho[x][y]=1;
      phi[x][y]=1;
      if (y-(ydim/10)<0) phi[x][y]=2;
      if (y-9*ydim/10>0) phi[x][y]=0.5;
      // phi[x][y]=1;
      u[x][y][0]=0;
      u[x][y][1]=0;
      if ((x==xdim/2)&&(y==ydim/2)) rho[x][y]=2;
    }
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++){
	  (*f)[i][j][x][y]=feq(i,j,rho[x][y],u[x][y][0],u[x][y][1]);
	  (*g)[i][j][x][y]=feq(i,j,phi[x][y],u[x][y][0],u[x][y][1]);
	}
  }
  defineBoundary();
}

void GetGraph(){
  if (rhosymreq){
    for (int x=0; x<xdim; x++)
      for (int y=0; y<ydim; y++)
	rhosym[x][y]=rho[x][y]-rho[y][x];
  }
  if (usymreq){
    for (int x=0; x<xdim; x++)
      for (int y=0; y<ydim; y++){
	usym[x][y][0]=u[x][y][0]-u[y][x][1];
	usym[x][y][1]=u[x][y][1]-u[y][x][0];
      }
  }
  if (fieldgreq){
    fieldgreq=0;
    for (int x=0; x<xdim; x++)
      for (int y=0; y<ydim; y++)
	fieldg[x][y]=field[x][y];
  }
}

double AveragePhiX(int y){
  double sum=0;
  for(int x; x<xdim; x++) sum+=phi[x][y];
  return sum/xdim;
}

void WriteAveragePhi(){

  FILE *res;
  char name[255];

  sprintf(name,"AvgPhi_UniformHeat_Tau_%e_Tauphiin_%e_Tauphiout_%e_Acc_%e_Heat_%e.dat",tau,tauphiin,tauphiout,Acc,heat);
  res = fopen(name,"w");
  for (int y=0; y<ydim; y++)
    fprintf(res,"%e %i %e\n",t,y,AveragePhiX(y));
  fclose(res);
}

int main(){
    struct timespec TimeSpec;

    int sstep=0,cont=0,done=0,Repeat=1;
  TimeSpec.tv_sec=0;
  TimeSpec.tv_nsec=10000000;
  DefineField();
  init();
  DefineGraphNxN_R("rho",&rho[0][0],&Xdim,&Ydim,NULL);
  DefineGraphNxN_R("phi",&phi[0][0],&Xdim,&Ydim,NULL);
  DefineGraphNxN_R("field",&fieldg[0][0],&Xdim,&Ydim,&fieldgreq);
  DefineGraphNxN_RxR("u",&u[0][0][0],&Xdim,&Ydim,NULL);
  DefineGraphNxN_R("rho sym",&rhosym[0][0],&Xdim,&Ydim,&rhosymreq);
  DefineGraphNxN_RxR("u sym",&usym[0][0][0],&Xdim,&Ydim,&usymreq);

  
  StartMenu("LB Diffusion",1);
  DefineGraph(contour2d_,"Graph");
  DefineFunction("Reinitialize",init);
  DefineFunction("Write Average Phi",WriteAveragePhi);
  DefineDouble("tau",&tau);
  DefineDouble("tauphiin",&tauphiin);
  DefineDouble("tauphiout",&tauphiout);
  DefineDouble("Acc",&Acc);
  DefineDouble("hole",&hole);
  DefineDouble("Heat",&heat);
  DefineDouble("Cool",&cool);
  DefineInt("Repeat",&Repeat);
  DefineBool("sstep",&sstep);
  DefineBool("cont",&cont);
  DefineBool("done",&done);
  EndMenu();

  while (!done){
    Events(1);
    GetGraph();
    DrawGraphs();
    nanosleep(&TimeSpec,NULL);
 
    if (sstep||cont){
      sstep=0;
      for (int i=0; i<Repeat; i++){
	iterate(*f,*ff,*g,*gg);
	farr *tmp; // switching the f's
	tmp=f;
	f=ff;
	ff=tmp;
 	tmp=g;
	g=gg;
	gg=tmp;
      }
    }
    
  }
}
