#include <stdio.h>
#include <stdlib.h>

int N=10;
double Nm=1000;
double a=0.1;
double b=0.65;

int D(int N){
  int Dtmp=0;
  double p=a*Nm/(Nm-N);
  if ((p<0)||(p>1)) p=1;
  for (int i=0; i<N; i++){
    if (rand()<p*RAND_MAX) Dtmp++;
  }
    return Dtmp;
}

int B(int N){
  int Btmp=0;
  for (int i=0; i<N; i++){
    if (rand()<b*RAND_MAX) Btmp++;
  }
    return Btmp;
}

int F(int N){
  int Bt=B(N);
  int Dt=D(N);
  return N+Bt-Dt;
}

int main(){
  FILE *res;

  res=fopen("ResDiscrete_b_0.65.dat","w");
  for (int i=0; i<100; i++){
    fprintf(res,"%i %i\n",i,N);
    for (int j=0; j<1; j++)
      N=F(N);
  }
}
