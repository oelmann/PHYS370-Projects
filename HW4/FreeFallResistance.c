#include <stdio.h>
#include <math.h>
#define D 1

double g=9.81;
double a=5;
double timelength=1;


void iterate(double x[D],double v[D],double dt){
  double F[D];
  for (int i=0; i<D; i++){
    F[i]=-g+a*v[i]*v[i];
    x[i]+=v[i]*dt;
    v[i]+=F[i]*dt;
  }
}

int main(){
  FILE *res;
  char name[100];
  double x0=0,v0=1;
  double x[D]={x0},v[D]={v0};
  double dt=0.01,t=0;
  
  sprintf(name,"FreeResistanceDt%f_a%f.dat",dt,a);
  res=fopen(name,"w");
  for (t=0; t<timelength; t+=dt){
    /* Energy = mgh+0.5 m v^2*/
    fprintf(res,"%f %f  %f %f\n",t,x[0],v[0],g*x[0]+0.5*v[0]*v[0]);
    iterate(x,v,dt);
  }
  fclose(res);

  sprintf(name,"FreeResistanceDt%f_a%fAnalytical.dat",dt,a);
  res=fopen(name,"w");
  for (t=0; t<timelength; t+=dt){
    fprintf(res,"%f %f  %f %f\n",t,x0+log(1+t*v0),v0/(1+t*v0),g*(x0+log(1+t*v0))+0.5*(v0/(1+t*v0))*(v0/(1+t*v0)));
  }
  fclose(res);

  sprintf(name,"TerminalVelocityDt%f_a%fAnalytical.dat",dt,a);
  res=fopen(name,"w");
  for (t=0; t<timelength; t+=dt){
    fprintf(res,"%f %f\n",t,-sqrt(g/a));
  }
  fclose(res);

}
