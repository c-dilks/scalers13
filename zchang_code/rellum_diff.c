#include<stdio.h>
#include<stdlib.h>

int main(int argc, char *argv[])
{
  if(argc!=3)
  {
    printf("the command line is *.exe *.dat *.dat\n");
    return 1;
  }

  FILE *fa,*fb;
  double ra[6],rb[6];
  double dr[6];

  if((fa=fopen(argv[1],"r"))==NULL||(fb=fopen(argv[2],"r"))==NULL)
  {
    printf("can't open %s or %s\n",argv[1],argv[2]);
    return 1;
  }
  int i;

  for(i=0;i<6;i++)
  {
    if((fscanf(fa,"%lf",&ra[i])==1)&&(fscanf(fb,"%lf",&rb[i])==1))
      dr[i]=ra[i]-rb[i];
    printf("%lf ",dr[i]);
  }

  printf("\n");

  return 0;
}
