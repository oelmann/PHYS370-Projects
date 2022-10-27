#include <stdio.h>
#include <math.h>

#define D 1

double k=10;
double totaltime = 10;
double x0=0,v0=1;

double dt=0.25,t=0;


void iterate1st(double x[D],double v[D],double dt){
  double F[D];
  for (int i=0; i<D; i++){
    F[i]=-k*x[i];
    x[i]+=v[i]*dt;
    v[i]+=F[i]*dt;
  }
}

void iterate2nd(double x[D],double v[D],double dt){
  double F[D];
  for (int i=0; i<D; i++){
    F[i]=-k*x[i];
    v[i]+=0.5*F[i]*dt;
    x[i]+=v[i]*dt;
    F[i]=-k*x[i];
    v[i]+=0.5*F[i]*dt;
  }
}

double sum2errorx(double x[D], double t, double om,  double error){
  double xth={x0*cos(om*t)+v0/om*sin(om*t)};
  double errorx = error + (x[0]-xth)*(x[0]-xth)/(totaltime/dt);
  return errorx;
}

double sum2errorv(double v[D], double t, double om, double error){
  double vth={-om*x0*sin(om*t)+v0*cos(om*t)};
  double errorv = error + (v[0]-vth)*(v[0]-vth)/(totaltime/dt);
  return errorv;
}

int main(){
  FILE *res;
  char name[100];
  double x[D]={x0},v[D]={v0};
  double errorx=0;
  //double errorv=0;
  double om=sqrt(k);

  
  sprintf(name,"HarmonicError%f_k%f_2nd.dat",dt,k);
  res=fopen(name,"w");
  for(double dt=0.0001; dt<0.5; dt=dt*1.1){
  x[0]=x0;v[0]=v0;
  errorx=0;//errorv=0;
  for (t=0; t<totaltime; t+=dt){
    errorx=sum2errorx(x,t,om,errorx);
    //errorv=sum2errorv(v,t,om,errorv);
    iterate2nd(x,v,dt);
  }
  fprintf(res,"%e %e ",dt,errorx);//,errorv);

  x[0]=x0;v[0]=v0;
  errorx=0;//errorv=0;
  

  for (t=0; t<totaltime; t+=dt){
    errorx=sum2errorx(x,t,om,errorx);
    //errorv=sum2errorv(v,t,om,errorv);
    iterate1st(x,v,dt);
  }
  fprintf(res,"%e \n",errorx);//,errorv);
  }
  
  fclose(res);
  


}
