#include <stdio.h>
#define D 1

double g=9.81;


void iterate(double x[D],double v[D],double dt){
  double F[D]={-g};
  for (int i=0; i<D; i++){
    x[i]+=v[i]*dt;
    v[i]+=F[i]*dt;
  }
}

int main(){
  FILE *res;
  char name[100];
  double x0=0,v0=1;
  double x[D]={x0},v[D]={v0};
  double dt=0.0025,t=0;
  
  sprintf(name,"FreeDt%f.dat",dt);
  res=fopen(name,"w");
  for (t=0; t<1; t+=dt){
    /* Energy = mgh+0.5 m v^2*/
    fprintf(res,"%f %f  %f %f\n",t,x[0],v[0],g*x[0]+0.5*v[0]*v[0]);
    iterate(x,v,dt);
  }
  fclose(res);

  sprintf(name,"FreeDt%fAnalytical.dat",dt);
  res=fopen(name,"w");
  for (t=0; t<1; t+=dt){
    fprintf(res,"%f %f  %f %f\n",t,x0+v0*t+0.5*g*t*t,v0+g*t,g*x0+0.5*v0*v0);
  }
  fclose(res);

}
