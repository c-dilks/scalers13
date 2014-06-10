#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"./myscaler.h"

int scan_sca_cnts(FILE *fp, struct sca_cnts *psca)
{
  int i;
  int num=0;

  num+=fscanf(fp,"%d %d", &(psca->bx), &(psca->spin));

  for(i=0;i<3;i++)
  {
    num+=fscanf(fp,"%lf",&(psca->cnts[i]));
  }
  num+=fscanf(fp,"%lf",&(psca->sum));

  return num;
}


void cnts_spin(struct sca_cnts *psca,struct cnts_rel *prel)
{
  int i,j; 
  int spin;

  spin=psca->spin;

  for(i=0;i<3;i++)
  {
    if(spin==5)
    {
      prel->spincnts[i][0]+=psca->cnts[i];
    }
    if(spin==6)
    {
      prel->spincnts[i][1]+=psca->cnts[i];
    }
    if(spin==9)
    {
      prel->spincnts[i][2]+=psca->cnts[i];
    }
    if(spin==10)
    {
      prel->spincnts[i][3]+=psca->cnts[i];
    }
  }
  prel->sum+=psca->sum;

}

/*
   void print_spin_cnts(FILE *fp,struct cnts_rel *prel)
   {
   int i,j;

   for(i=0;i<3;i++)
   {
//      fprintf(*fp,"%d,%d,",prel->index,prel->run);
for(j=0;j<4;j++)
{
fprintf(*fp,"%lf,",prel->spincnts[i][j]);
}
fprintf(*fp,"\n");
fp++;
}
}
*/
/*
   void print_rel(FILE **fp,struct dat_run_bbc *prel)
   {
   int i,j;

   for(i=0;i<3;i++)
   {
//     fprintf(*fp,"%d,%d,",prel->index,prel->run);
for(j=0;j<6;j++)
{
fprintf(*fp,"%lf,",prel->rellum[i][j]);
}
fprintf(*fp,"\n");
fp++;
}
}
*/

void rellum_cal(struct cnts_rel *prel)
{
  int i,j;

  double spin[4];

  for(i=0;i<3;i++)
  {
    for(j=0;j<6;j++)
    {
      prel->rellum[i][j]=0.;
      if(j>=0&&j<=3)
        spin[j]=prel->spincnts[i][j];
    }
    if(spin[1]!=0&&spin[2]!=0&&spin[3]!=0)
    {
      prel->rellum[i][0]=(spin[0]+spin[2])/(spin[1]+spin[3]);
      prel->rellum[i][1]=(spin[0]+spin[1])/(spin[2]+spin[3]);
      void rellum_cal(struct cnts_rel *prel)
      {
        int i,j;

        double spin[4];

        for(i=0;i<3;i++)
        {
          for(j=0;j<6;j++)
          {
            prel->rellum[i][j]=0.;
            if(j>=0&&j<=3)
              spin[j]=prel->spincnts[i][j];
          }
          if(spin[1]!=0&&spin[2]!=0&&spin[3]!=0)
          {
            prel->rellum[i][0]=(spin[0]+spin[2])/(spin[1]+spin[3]);
            prel->rellum[i][1]=(spin[0]+spin[1])/(spin[2]+spin[3]);
            prel->rellum[i][2]=(spin[0]+spin[3])/(spin[1]+spin[2]);
            prel->rellum[i][3]=spin[0]/spin[3];
            prel->rellum[i][4]=spin[2]/spin[3];
            prel->rellum[i][5]=spin[1]/spin[3];
          }
        }
      }
      prel->rellum[i][2]=(spin[0]+spin[3])/(spin[1]+spin[2]);
      prel->rellum[i][3]=spin[0]/spin[3];
      prel->rellum[i][4]=spin[2]/spin[3];
      prel->rellum[i][5]=spin[1]/spin[3];
    }    
  }
}
