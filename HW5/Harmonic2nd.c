#include <stdio.h>
#include <math.h>

#define D 1

double k=10;


void iterate(double x[D],double v[D],double dt){
  double F[D];
  for (int i=0; i<D; i++){
    F[i]=-k*x[i];
    v[i]+=0.5*F[i]*dt;
    x[i]+=v[i]*dt;
    F[i]=-k*x[i];
    v[i]+=0.5*F[i]*dt;
  }
}

int main(){
  FILE *res;
  char name[100];
  double x0=0,v0=1;
  double x[D]={x0},v[D]={v0};
  double dt=0.1,t=0;

  double om=sqrt(k);
  
  sprintf(name,"HarmonicDt%f_2nd.dat",dt);
  res=fopen(name,"w");
  for (t=0; t<1; t+=dt){
    /* Energy = mgh+0.5 m v^2*/
    fprintf(res,"%f %f  %f %f\n",t,x[0],v[0],0.5*k*x[0]*x[0]+0.5*v[0]*v[0]);
    iterate(x,v,dt);
  }
  fclose(res);

  sprintf(name,"HarmonicDt%fAnalytical2nd.dat",dt);
  res=fopen(name,"w");
  for (t=0; t<1; t+=dt){
    fprintf(res,"%f %f  %f %f\n",t,x0*cos(om*t)+v0/om*sin(om*t),-om*x0*sin(om*t)+v0*cos(om*t),0.5*k*x0*x0+0.5*v0*v0);
  }
  fclose(res);

}
