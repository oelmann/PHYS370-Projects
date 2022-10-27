#include<stdio.h>
#include<math.h>

#define Count 1000000

double Ninit = 1;
double Nm = 100;
double a = 0.1;
int count=0;
double dat[Count][3];

double binit=0.4, bend=0.65, bstep=0.001;
int Nequalize=1000, Nshow=100;

double F(double N,double b){
  return N+b*N-a*N*Nm/(Nm-N);
}

double Fderive(double x0, double b){
  const double delta = 1;
  double x1 = x0 - delta;
  double x2 = x0 + delta;
  double y1 = F(x1,b);
  double y2 = F(x2,b);
  return (y2-y1) / (x2-x1);  
}

void GetData(){
   count=0;
   for (double b=binit; b<bend; b+=bstep){
    double N=Ninit;
    for (int i=0; i<Nequalize; i++){
      N = F(N,b);
      
  }
    for (int i=0; i<Nshow; i++){
      N = F(N,b);
      if (N<0) N=0;
      dat[count][0]=b;
      dat[count][1]=N;
      dat[count][2]=Fderive(dat[count][1],b);
      count = count + 1;
      if (count==Count){
	count = count - 1;
	printf("Increase value of Count to allow more points!\n");
      }
    }
  }
}

void WriteToFile(){

  FILE *res;
  char name[255];

  sprintf(name,"BunnyLatebinit_%e_bstop_%e.dat",binit,bend);
  res = fopen(name,"w");
  for (int i=0; i<count; i++)
    fprintf(res,"%e %e\n",dat[i][0],dat[i][1]);
  fclose(res);
}

  int main(){
    GetData();
    WriteToFile();

  FILE *res;
  
  res = fopen("N_vs_Fderiv.dat","w");
  for (int i=0; i<count; i++){
    fprintf(res,"%e %e\n",dat[i][1],dat[i][2]);
  fclose(res);
  }
 
}
