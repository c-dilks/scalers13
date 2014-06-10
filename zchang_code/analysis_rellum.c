#include<stdio.h>
#include<stdlib.h>

#include"./myscaler.h"

int main(int argc, char *argv[])
{
  if(argc!=2)
  {
    printf("the command line is *.exe *.dat\n");
    return 1;
  }

  FILE *inf;
  FILE *outfrel[3];
  FILE *outfcnts[3];
  int i,j;

  // structs defined in myscaler.h
  struct sca_cnts scaCnts;
  struct cnts_rel cntsRel;

  char tag[4];
  char filerel[3][50];
  char filecnts[3][50];

  char *pos[3]={"e","w","x"};
  //  int run;

  // open datfile
  if((inf=fopen(argv[1],"r"))==NULL)
  {
    printf("can't open %s\n",argv[1]);
    return 1;
  }

  //  sscanf(argv[1],"%s_run%d_4.dat",&tag,&run);

  // set up output files for writing
  for(i=0;i<3;i++)
  {
    sprintf(filerel[i],"rel_%s_%s",pos[i],argv[1]);
    outfrel[i]=fopen(filerel[i],"w");

    sprintf(filecnts[i],"cnts_%s_%s",pos[i],argv[1]);
    outfcnts[i]=fopen(filecnts[i],"w");
  }

  // zero cntsRel struct spincnts entries
  for(i=0;i<3;i++)
  {
    for(j=0;j<4;j++)
    {
      cntsRel.spincnts[i][j]=0.;
    }
  }

  // scans datfile (six columns)
  // - bx
  // - spin
  // - cnts[0,1,2]
  // - sum
  while(scan_sca_cnts(inf,&scaCnts)==6)
  {
    cnts_spin(&scaCnts,&cntsRel);
  }
  fclose(inf);


  rellum_cal(&cntsRel);

  for(i=0;i<3;i++)
  {
    for(j=0;j<6;j++)
    {
      printf("%lf ", cntsRel.rellum[i][j]);
      fprintf(outfrel[i],"%lf ",cntsRel.rellum[i][j]);
    }
    printf("\n");
    fprintf(outfrel[i], "\n");
    fclose(outfrel[i]);

    for(j=0;j<4;j++)
    {
      //	  printf("%lf ",cntsRel.spincnts[i][j]);	      
      fprintf(outfcnts[i],"%lf ",cntsRel.spincnts[i][j]);
    }

    //      printf("\n");
    fprintf(outfcnts[i],"\n");
    fclose(outfcnts[i]);
  }

  return 0;

}
