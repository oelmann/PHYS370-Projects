#include<stdio.h>
#include<mygraph.h>

#define Count 1000000

double Ninit = 1;
double Nm = 100;
double a = 0.1;
int count=0;
double dat[Count][2];
int points = 100;
double sol1[100][2];
double sol2[100][2];
double temp_sol2[100][2];
double breakpoint[2];
double sol3[100][2];

double binit=0.4, bend=0.65, bstep=0.001;
int Nequalize=1000, Nshow=100;

//This is function that models new population based on current population N and birth rate b
double F(double N,double b){
  return N+b*N-a*N*Nm/(Nm-N);
}

//This function solves for the fixed point values
void GetSol1(){
  for (int i=0; i<100; i++){
    double b=binit + (bend - binit)*i/100;
    sol1[i][0]=b;
    sol1[i][1]=(b-a)*Nm/b;
  }
}

//This function solves for the slope of the fixed points and finds where bifurcation occurs
void GetSol2(){
  for (int i=0; i<100; i++){
    double b=binit + (bend - binit)*i/100;
    double x = (b-a)*Nm/b;
    double alpha = Nm/(Nm-x);
    sol2[i][0]=b;
    sol2[i][1]=(1+b-a*alpha)-x*a*alpha/(Nm-x);
  }
  int i=0;
  while((sol2[i][1] < 1) && (sol2[i][1] > -1)){
    temp_sol2[i][0] = sol2[i][0];
    temp_sol2[i][1] = sol2[i][1];
    i++;
  }
  GetSol1();
  for (int j=0; j<100; j++){
    sol2[j][0] = sol2[i][0];
    sol2[j][1] = sol1[i][1];
    //printf("%e %e\n",sol2[j][0], sol2[j][1]);
  }
}


void GetSol3(){
  for (int i=0; i<100; i++){
    double b=binit + (bend - binit)*i/100;
    sol3[i][0]=b;
    sol3[i][1]=(-a*Nm+b*Nm+Nm)/(b+1);
  }
}


void GetData(){
   count=0;
   for (double b=binit; b<bend; b+=bstep){
    double N=Ninit;
    //This loop takes the population to its max for a given b
    for (int i=0; i<Nequalize; i++){
      N = F(N,b);
  }
    for (int i=0; i<Nshow; i++){
      N = F(N,b);
      if (N<0) N=0;
      dat[count][0]=b;
      dat[count][1]=N;
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
    int done=0;


    SetDefaultLineType(0);
    SetDefaultShape(4);
    DefineGraphN_RxR("BD",&dat[0][0],&count,NULL);
    SetDefaultColor(2);
    SetDefaultLineType(1);
    SetDefaultShape(0);
    DefineGraphN_RxR("Sol1",&sol1[0][0],&points,NULL);
    SetDefaultShape(3);
    SetDefaultLineType(0);
    DefineGraphN_RxR("Sol2",&sol2[0][0],&points,NULL);
    SetDefaultLineType(1);
    SetDefaultShape(0);
    DefineGraphN_RxR("Sol3",&sol3[0][0],&points,NULL);
  
    
    StartMenu("Bunny Graphs",1);
    DefineDouble("a",&a);
    DefineDouble("binit",&binit);
    DefineDouble("bend",&bend);
    DefineDouble("bstep",&bstep);
    DefineInt("Nequalize",&Nequalize);
    DefineInt("Nshow",&Nshow);
    DefineFunction("GetData",GetData);
    DefineFunction("GetSol1",GetSol1);
    DefineFunction("GetSol2",GetSol2);
    DefineFunction("GetSol3",GetSol3);
    DefineGraph(curve2d_,"Data Graph");
    DefineBool("done",&done);
    EndMenu();

    while(done==0){
      Events(1);
      DrawGraphs(1);
    }
 
}
