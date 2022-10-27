#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mygraph.h>

#define D 3
#define N 27

double Dim[D];
double x0[N][D],v0[N][D],M[N];
double x[N][D],v[N][D];
double dt=0.01,t=0;
// GUI variables
double xcenter=0,ycenter=0,Scale=600;
double screend=1,objd=10, R=0.3,tet=0;
typedef struct part_ {double x[D]; int c;} part;

double T=0,Tmeas;
#define L 10000 
double trace[N][L][D];
int traceL=0,Lmax=L;

void draw3dtrace(int xdim, int ydim){
  double tracel[N][L][D];
  /*  for (int n=0; n<N; n++){
    xl[n].c=n+1;
    for (int d=0; d<D; d++)
      xl[n].x[d]=x[n][d];
      }*/ 
   memcpy(&tracel[0][0][0],&trace[0][0][0],N*L*D*sizeof(double));
  // rotate around y axis
  double ct=cos(tet);
  double st=sin(tet);
  for (int n=0; n<N; n++){
    for (int l=0; l<traceL; l++){
      double ztmp= ct*tracel[n][l][2]+st*tracel[n][l][0];
      double xtmp=-st*tracel[n][l][2]+ct*tracel[n][l][0];
      tracel[n][l][0]=xtmp;
      tracel[n][l][2]=ztmp;
    }
  }
  
  // qsort(&xl[0].x[0],N,sizeof(part),compare);
  for (int n=0; n<N; n++){
    XPoint xp[traceL];
    int l=0;

    while (l<traceL){
      int broken=0;
      int nop=0;
      xp[nop].x=xdim/2+Scale*screend*tracel[n][l][0]/(tracel[n][l][2]+objd);
      xp[nop].y=ydim/2-Scale*screend*tracel[n][l][1]/(tracel[n][l][2]+objd);
      l++;
      nop++;
      while ((broken==0)&&(l<traceL)){
	xp[nop].x=xdim/2+Scale*screend*tracel[n][l][0]/(tracel[n][l][2]+objd);
	xp[nop].y=ydim/2-Scale*screend*tracel[n][l][1]/(tracel[n][l][2]+objd);
	for (int d=0; d<D; d++){
	  if(fabs(tracel[n][l][d]-tracel[n][l-1][d])>Dim[d]/2) broken=1;
	}
	if (broken==0){
	  l++;
	  nop++;
	}
      }
      int col=0;
      if (n<N/2) 
	col=2;
      else col=3;
      
      myline_polygon(col,xp,nop);
    }
  }
}
 
int compare(part *x1,part *x2){
  if (x1->x[2]<x2->x[2]) return 1;
  else return -1;
}

void draw3d(int xdim, int ydim){
  part xl[N];
  for (int n=0; n<N; n++){
    if (n<N/2) 
      xl[n].c=2;
    else xl[n].c=3;
    for (int d=0; d<D; d++)
      xl[n].x[d]=x[n][d]-Dim[d]*0.5;
  }
  //  memcpy(&xl[0].x[0],&x[0][0],N*sizeof(part));
  // rotate around y axis

  double ct=cos(tet);
  double st=sin(tet);
  for (int n=0; n<N; n++){
    double ztmp= ct*xl[n].x[2]+st*xl[n].x[0];
    double xtmp=-st*xl[n].x[2]+ct*xl[n].x[0];
    xl[n].x[0]=xtmp;
    xl[n].x[2]=ztmp;
  }
  
  qsort(&xl[0].x[0],N,sizeof(part),compare);
  for (int n=0; n<N; n++){
  int xc=xdim/2+Scale*screend*xl[n].x[0]/(xl[n].x[2]+objd);
    int yc=ydim/2-Scale*screend*xl[n].x[1]/(xl[n].x[2]+objd);
    int rc=Scale*R*screend/(xl[n].x[2]+objd);
    myfilledcircle(xl[n].c,xc,yc,rc);
  }
}

void draw(int xdim, int ydim){
  for (int n=0; n<N; n++){
    int xc=xdim/2+(x[n][0]-xcenter)*Scale;
    int yc=ydim/2-(x[n][1]-ycenter)*Scale;
    int rc=10;
    myfilledcircle(n+1,xc,yc,rc);
  }
}
void draw_trace(int xdim, int ydim){
  if (traceL>1){
    for (int n=0; n<N; n++){
      XPoint xp[traceL];
      for (int l=0; l<traceL; l++){
	xp[l].x=xdim/2+(trace[n][l][0]-xcenter)*Scale;
	xp[l].y=ydim/2-(trace[n][l][1]-ycenter)*Scale;
      }
      myline_polygon(n+1,xp,traceL);
    }
  }
}

void update_trace(){
  while (traceL>=Lmax){
    traceL--;
    for (int n=0; n<N; n++){
      for (int l=0; l<traceL; l++)
	for (int d=0; d<D; d++)
	  trace[n][l][d]=trace[n][l+1][d];
    }
  }
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      trace[n][traceL][d]=x[n][d]-Dim[d]/2;
    }
  }
  traceL++;
}

void getP(double v[N][D],double M[N],double P[D]){
  for (int d=0; d<D; d++) P[d]=0;
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      P[d]+=M[n]*v[n][d];
    }
  }
}

void SetP0(double v[N][D],double M[N]){
  double P[D];
  getP(v,M,P);
  double MM=0;
  for (int n=0; n<N; n++) MM+=M[n];
  double V[D];
  for (int d=0; d<D; d++) V[d]=P[d]/MM;

  for (int n=0; n<N; n++)
    for (int d=0; d<D; d++)
      v[n][d]-=V[d];
}

double getE_P(double M[N],double x[N][D],double v[N][D],double P[D]){
  // Ek+Ep
  double E=0;
  for (int d=0; d<D; d++) P[d]=0;
  for (int n=0; n<N; n++){
    double v2=0;
    for (int d=0; d<D; d++){
      P[d]+=M[n]*v[n][d];
      v2+=v[n][d]*v[n][d];
    }
    E+=0.5*M[n]*v2;
  }

  // Ep = -G*m1*m2/r
  for (int n=0; n<N; n++){
    for (int m=0; m<n; m++){
      double dr[D][3];
      for (int d=0; d<D; d++){
	dr[d][1]=x[m][d]-x[n][d];
	dr[d][2]=dr[d][1]+Dim[d];
	dr[d][0]=dr[d][1]-Dim[d];
      }
      for (int k=0; k<pow(3,D); k++){ // funky way of summing over all neighbors in D dimensions.
	int kk=k;
	int ind[D];
	for (int d=0; d<D; d++){
	  ind[d]=kk%3;
	  kk/=3;
	}
	double r2=0;
	for (int d=0; d<D; d++){
	  r2+=dr[d][ind[d]]*dr[d][ind[d]];
	}
	double r6=r2*r2*r2;
	E+=1/(r6*r6)-1/(r6);
      }
    }
  }
  return E;
}


void getF(double F[N][D]){
  double dr[D][3],absdr2;
  memset(&F[0][0],0,N*D*sizeof(double));
  // this should be for arbitrary numbers of dimensions.. but this is only for 3 dimensions
  for (int n=0; n<N; n++)
    for (int m=0; m<n; m++){
      for (int d=0; d<D; d++){
	dr[d][1]=x[m][d]-x[n][d];
	dr[d][2]=dr[d][1]+Dim[d];
	dr[d][0]=dr[d][1]-Dim[d];
      }
      
      for (int k=0; k<pow(3,D); k++){ // funky way of summing over all neighbors in D dimensions.
	int kk=k;
	int ind[D];
	for (int d=0; d<D; d++){
	  ind[d]=kk%3;
	  kk/=3;
	}
	absdr2=0;
	for (int d=0; d<D; d++){
	  absdr2+=dr[d][ind[d]]*dr[d][ind[d]];
	}
	for (int d=0; d<D; d++){
	  // Lennard Jones potential 1/r^12-1/r^6
	  // Force = -dV/dr = (12/r^13 - 6/r^7) vec(r)/r
	  double r4=absdr2*absdr2;
	  double Fscal=12/(r4*r4*r4*absdr2)-6/(r4*r4);
	  F[n][d]-=Fscal*dr[d][ind[d]];
	  F[m][d]+=Fscal*dr[d][ind[d]];
	}
      }
    }
}

void iterate(double x[N][D],double v[N][D],double dt){
  double F[N][D]; // \vec{F} = G m*m/r^3 \vec{r}
  // get the force on all N particles
  getF(F);

  // update positions and velocities with the Verlet algorithm
  for (int n=0; n<N; n++){
    for (int i=0; i<D; i++){
      v[n][i]+=0.5*F[n][i]/M[n]*dt;
      x[n][i]+=v[n][i]*dt;
      x[n][i]=fmod(x[n][i]+Dim[i],Dim[i]); // Periodic boundary [0,Dim[i][
    }
  }
  getF(F);
  for (int n=0; n<N; n++){
    for (int i=0; i<D; i++){
      v[n][i]+=0.5*F[n][i]/M[n]*dt;
    }
  }
}

double GetT(){
  // this assumes that the mean velocity is zero
  double kbT=0;
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      kbT+=M[n]*v[n][d]*v[n][d];
    }
  }
  return kbT/(N*D);
}



void SetPzero(){
  double p[D],m=0;
  for (int d=0; d<D; d++) p[d]=0;
  
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      p[d]+=M[n]*v[n][d];
    }
    m+=M[n];
  }
  printf("P = (%e; %e; %e)\n ",p[0],p[1],p[2]);
  for (int d=0; d<D; d++) p[d]/=m;
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      v[n][d]-=p[d];
    }
  }
  for (int d=0; d<D; d++) p[d]=0;
  
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      p[d]+=M[n]*v[n][d];
    }
  }
  printf("P = (%e; %e; %e)\n ",p[0],p[1],p[2]);
}

void SetT(){
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      v[n][d]= 2*((double)rand()/RAND_MAX-0.5)*sqrt(3*T/M[n]);
    }
  }
  SetPzero();
}

void initconfig(){
  for (int d=0; d<D; d++) Dim[d]=10;
  int n=0;
  int LN=3;
  for (int x=0; (x<LN)&&(n<N); x++)
    for (int y=0; (y<LN)&&(n<N); y++)
      for (int z=0; (z<LN)&&(n<N); z++){
	x0[n][0]=x-0.5+Dim[0]/2;
	x0[n][1]=y-0.5+Dim[1]/2;
	x0[n][2]=z-0.5+Dim[2]/2;
	n++;
      }
  for (int n=0; n<N; n++){
    if (n<N/2) M[n]=1;
    else M[n]=5;
    for (int d=0; d<D; d++){
      v0[n][d]=0; //(double)rand()/RAND_MAX;
    }
  }
}

void WriteV(double v[N][D]){
  FILE *res = fopen("MeasuredVelocity.dat","w");

  
  for (int n=0; n<N; n++)
    fprintf(res,"%i %e\n",i,);
  fclose(res);
} 

void init(){
  t=0;
  traceL=0;
  for (int n=0; n<N; n++)
    for (int i=0; i<D; i++){
      x[n][i]=x0[n][i];
      v[n][i]=v0[n][i];
    }
  SetP0(v,M);
}

int main(){
  struct timespec TimeSpec;
  int done=0, cont=0,sstep=0,repeat=1;
  double E,P[D];
  AddFreedraw("planets",draw);
  AddFreedraw("traces",draw_trace);
  if (D==3){
    AddFreedraw("3d trace",draw3dtrace);
    AddFreedraw("3d view",draw3d);  
  }
  StartMenu("Planets",1);
  StartMenu("Init",0);
  /*
  for (int n=0; n<N; n++){
    DefineDouble("m",&M[n]);
    for (int d=0; d<D; d++){
      DefineDouble("x",&x0[n][d]);
    }
    for (int d=0; d<D; d++){
      DefineDouble("v",&v0[n][d]);
    }
  }
  */
  DefineFunction("Init",init);
  EndMenu();
  DefineGraph(freedraw_,"Planets");
  DefineDouble("Scale",&Scale);
  if (D==3){
    DefineDouble("screen dist",&screend);
    DefineDouble("obj dist",&objd);
    DefineDouble("R",&R);
    DefineDouble("tety",&tet);
  }
  DefineMod("Trace length",&Lmax,L);
  DefineLong("wait",&TimeSpec.tv_nsec);
  DefineDouble("dt",&dt);
  DefineDouble("Tmeas",&Tmeas);
  DefineDouble("T",&T);
  DefineFunction("SetT",SetT);
  DefineDouble("E",&E);
  for (int d=0; d<D; d++) DefineDouble("P",&P[d]);
  DefineInt("Repeat",&repeat);
  DefineBool("sstep",&sstep);
  DefineBool("cont",&cont);
  DefineBool("done",&done);
  EndMenu();

  initconfig();
  init();
  TimeSpec.tv_sec=0;
  TimeSpec.tv_nsec=10000000;
  
  while (!done){
    Tmeas=GetT();
    Events(1);
    DrawGraphs();
    nanosleep(&TimeSpec,NULL);
    E=getE_P(M,x,v,P);
    if (cont||sstep){
      sstep=0;
      update_trace();
      for (int i=0; i<repeat; i++){
	iterate(x,v,dt);
	t+=dt;
      }
    }
  }
}
