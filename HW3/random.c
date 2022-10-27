#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N = 100;
int cycles = 1000;

double Fact(int N){
  double res=1;
  for(int i=2; i<=N; i++) res *=i;
  return res;
}

double Binomial(int N, double p){
    int success = 0;
    for(int i=0; i<N; i++)
      if(rand()<RAND_MAX*p) success++;
    return success;
}

double IdealBinomial(int N, double p, int i){
  double success = Fact(N)/Fact(i)/Fact(N-i)*pow(p,i)*pow(1-p,N-i);
  return success;
}


void WriteBinomial1(){
    FILE *res = fopen("binomial1.dat","w");
  for (int i=0; i<=cycles; i++)
    fprintf(res,"%e\n",Binomial(100,0.5));
  fclose(res);
}

void WriteBinomial2(){
    FILE *res = fopen("binomial2.dat","w");
  for (int i=0; i<=cycles; i++)
    fprintf(res,"%e\n",Binomial(N,0.2));
  fclose(res);
}

void WriteIdealBinomial1(){
  FILE *res = fopen("idealbinomial1.dat","w");
  for (int i=0; i<=100; i++)
    fprintf(res,"%i %e\n",i,IdealBinomial(100,0.5,i));
  fclose(res);
} 

void WriteIdealBinomial2(){
  FILE *res = fopen("idealbinomial2.dat","w");
  for (int i=0; i<=100; i++)
    fprintf(res,"%i %e\n",i,IdealBinomial(100,0.2,i));
  fclose(res);
} 

int main(){
    WriteBinomial1();
    WriteBinomial2();
    WriteIdealBinomial1();
    WriteIdealBinomial2();
}        
