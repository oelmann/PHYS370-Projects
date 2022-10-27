#include<stdio.h>
#include<math.h>

int main(){
  FILE *res;

  res = fopen("sin.dat","w");
  for (int i=0; i<100; i++){
    double x = i / (2*3.14159);
    fprintf(res,"%e %e\n",x,sin(x));
    x=x/i;
  }
  fclose(res);
  return 0;
}
