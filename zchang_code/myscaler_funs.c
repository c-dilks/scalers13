#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"myscaler.h"

int scan_sca_val(FILE *fp,struct sca_val *psca)
{
  int i,num=0;

  num+=fscanf(fp,"%d",&(psca->bx));

  for(i=0;i<8;i++)
  {
    num+=fscanf(fp,"%lld",&(psca->valBBC[i]));
  }
  for(i=0;i<8;i++)
  {
    num+=fscanf(fp,"%lld",&(psca->valZDC[i]));
  }
  for(i=0;i<4;i++)
  {
    num+=fscanf(fp,"%lld",&(psca->valVPD[i]));
  }
  num+=fscanf(fp,"%lld",&(psca->valSum));

  return num;
}
/*
   int scan_dat_bbc(FILE *fp,struct dat_bbc *psca)
   {
   int i,num=0;

   num+=fscanf(fp,"%d,%d,",&(psca->bx),&(psca->spin));

   for(i=0;i<3;i++)
   {
   num+=fscanf(fp,"%lf,",&(psca->nbbc[i]));
   }
   num+=fscanf(fp,"%lf",&(psca->sum));

   return num;
   }
   */

void print_sca_cnts(FILE *fp,struct sca_cnts *dat_appr)
{
  int i;

  fprintf(fp,"%d %d ",dat_appr->bx,dat_appr->spin);

  for(i=0;i<3;i++)
  {
    fprintf(fp,"%lf ",dat_appr->cnts[i]);
  }
  fprintf(fp,"%lf\n",dat_appr->sum);
}


void cnts_bbc(struct sca_val *pval,struct sca_cnts *pcnts)
{
  pcnts->cnts[0]=pval->valBBC[1]+pval->valBBC[3]+pval->valBBC[5]+pval->valBBC[7]; // east
  pcnts->cnts[1]=pval->valBBC[2]+pval->valBBC[3]+pval->valBBC[6]+pval->valBBC[7]; // west
  pcnts->cnts[2]=pval->valBBC[3]+pval->valBBC[7]; // coincidence

  pcnts->sum=pval->valSum;
}
void cnts_zdc(struct sca_val *pval,struct sca_cnts *pcnts)
{
  pcnts->cnts[0]=pval->valZDC[1]+pval->valZDC[3]+pval->valZDC[5]+pval->valZDC[7]; // east
  pcnts->cnts[1]=pval->valZDC[2]+pval->valZDC[3]+pval->valZDC[6]+pval->valZDC[7]; // west
  pcnts->cnts[2]=pval->valZDC[3]+pval->valZDC[7]; // coincidence

  pcnts->sum=pval->valSum;
}
void cnts_vpd(struct sca_val *pval,struct sca_cnts *pcnts)
{
  pcnts->cnts[0]=pval->valVPD[1]+pval->valVPD[3]; // east
  pcnts->cnts[1]=pval->valVPD[2]+pval->valVPD[3]; // west
  pcnts->cnts[2]=pval->valVPD[3]; // coincidence

  pcnts->sum=pval->valSum;
}
void assign_spin(FILE *fp, int *spinbit)
{
  int bx,spin;

  while(fscanf(fp,"%d %d",&bx,&spin)==2)
  {
    spinbit[bx]=spin;
  }
}
/*
   void cnts_appr_bbc(char *appr,struct sca_bx *pbx,struct dat_bbc *pbbc)
   {
   double bbce,bbcw,bbcx;
   double sum;

   sum=pbx->val;

   if(strcmp(appr,"na")==0)
   {
   printf("Calculating A/M Corrected Counts by na approach\n");

   bbce=pbx->valbbc[1]+pbx->valbbc[7];
   bbcw=pbx->valbbc[2]+pbx->valbbc[7];
   bbcx=pbx->valbbc[7];
   }

   if(strcmp(appr,"ga")==0)
   {
   printf("Calculating A/M Corrected Counts by ga approach\n");

   bbce=pbx->valbbc[1]+pbx->valbbc[3]+pbx->valbbc[4]+pbx->valbbc[5]+pbx->valbbc[6]+pbx->valbbc[7];
   bbcw=pbx->valbbc[2]+pbx->valbbc[3]+pbx->valbbc[4]+pbx->valbbc[5]+pbx->valbbc[6]+pbx->valbbc[7];
   bbcx=pbx->valbbc[4]+pbx->valbbc[5]+pbx->valbbc[6]+pbx->valbbc[7];
   }

   if(strcmp(appr,"sga")==0)
   {
   printf("Calculating A/M Corrected Counts by sga approach\n");

   bbce=pbx->valbbc[1]+pbx->valbbc[3]+pbx->valbbc[4]+pbx->valbbc[5]+pbx->valbbc[6]+pbx->valbbc[7];
   bbcw=pbx->valbbc[2]+pbx->valbbc[3]+pbx->valbbc[4]+pbx->valbbc[5]+pbx->valbbc[6]+pbx->valbbc[7];
   bbcx=pbx->valbbc[3]+pbx->valbbc[4]+pbx->valbbc[5]+pbx->valbbc[6]+pbx->valbbc[7];

   }

   if(strcmp(appr,"tda")==0)
   {
   printf("Calculating A/M Corrected Counts by tda approach\n");

   bbce=pbx->valbbc[1]+pbx->valbbc[3]+pbx->valbbc[7];
   bbcw=pbx->valbbc[2]+pbx->valbbc[3]+pbx->valbbc[7];
   bbcx=pbx->valbbc[3]+pbx->valbbc[7];
   }

   if(strcmp(appr,"mtda")==0)
   {   
   printf("Calculating A/M Corrected Counts by mtda approach\n");

   bbce=pbx->valbbc[1]+pbx->valbbc[3]+pbx->valbbc[5]+pbx->valbbc[7];
   bbcw=pbx->valbbc[2]+pbx->valbbc[3]+pbx->valbbc[6]+pbx->valbbc[7];
   bbcx=pbx->valbbc[3]+pbx->valbbc[7];
   }


   pbbc->nbbc[0]=bbce;
   pbbc->nbbc[1]=bbcw;
   pbbc->nbbc[2]=bbcx;

   pbbc->sum=sum;

   }
   */
void accidental(struct sca_cnts *dat_appr, struct sca_cnts *dat_base)
{
  double bbce,bbcw,bbcx;
  double sum;

  dat_appr->bx=dat_base->bx;
  dat_appr->spin=dat_base->spin;

  sum=dat_base->sum;

  bbce=dat_base->cnts[0];
  bbcw=dat_base->cnts[1];
  bbcx=dat_base->cnts[2];

  dat_appr->cnts[0]=(bbce-bbcx)/(sum-bbcw)*sum;
  dat_appr->cnts[1]=(bbcw-bbcx)/(sum-bbce)*sum;
  dat_appr->cnts[2]=(bbcx-bbce*bbcw/sum)/(sum+bbcx-bbce-bbcw)*sum;

  dat_appr->sum=sum;  
}

void multiple(struct sca_cnts *dat_mult,struct sca_cnts *dat_base)
{
  double bbce,bbcw,bbcx;
  double sum;

  sum=dat_base->sum;

  bbce=dat_base->cnts[0];
  bbcw=dat_base->cnts[1];
  bbcx=dat_base->cnts[2];

  dat_mult->bx=dat_base->bx;
  dat_mult->spin=dat_base->spin;

  dat_mult->cnts[0]=(-1)*sum*log(1.-bbce/sum);
  dat_mult->cnts[1]=(-1)*sum*log(1.-bbcw/sum);
  dat_mult->cnts[2]=(-1)*sum*log(1.-bbcx/sum);

  dat_mult->sum=sum;
}
/*
   void cnts_spin_bbc(FILE *fp,struct dat_run_bbc *pbbcr)
   {
   int i,j; 
   int spin;

//  void check_fill(struct dat_bbc *dap,int fill);

for(i=0;i<3;i++)
{
for(j=0;j<4;j++)
{
pbbcr->spincnts[i][j]=0.;
}
}
pbbcr->sum=0.;
while(scan_dat_bbc(fp,&(pbbcr->datbbc))==6)
{
spin=pbbcr->datbbc.spin;

//  check_fill(&(pbbcr->datbbc),pbbcr->fill);

for(i=0;i<3;i++)
{
if(spin==5)
{
pbbcr->spincnts[i][0]+=pbbcr->datbbc.cnts[i];
}
if(spin==6)
{
pbbcr->spincnts[i][1]+=pbbcr->datbbc.cnts[i];
}
if(spin==9)
{
pbbcr->spincnts[i][2]+=pbbcr->datbbc.cnts[i];
}
if(spin==10)
{
pbbcr->spincnts[i][3]+=pbbcr->datbbc.cnts[i];
}

}
pbbcr->sum+=pbbcr->datbbc.sum;
}
}

void print_spin_cnts(FILE **fp,struct dat_run_bbc *pbbcr)
{
int i,j;

for(i=0;i<3;i++)
{
fprintf(*fp,"%d,%d,",pbbcr->index,pbbcr->run);
for(j=0;j<4;j++)
{
fprintf(*fp,"%lf,",pbbcr->spincnts[i][j]);
}
fprintf(*fp,"\n");
fp++;
}
}

void print_rel(FILE **fp,struct dat_run_bbc *pbbcr)
{
int i,j;

for(i=0;i<3;i++)
{
fprintf(*fp,"%d,%d,",pbbcr->index,pbbcr->run);
for(j=0;j<6;j++)
{
fprintf(*fp,"%lf,",pbbcr->rellum[i][j]);
}
fprintf(*fp,"\n");
fp++;
}
}

void rellum_cal(struct dat_run_bbc *pbbcr)
{
  int i,j;

  double spin[4];

  for(i=0;i<3;i++)
  {
    for(j=0;j<6;j++)
    {
      pbbcr->rellum[i][j]=0.;
      if(j>=0&&j<=3)
        spin[j]=pbbcr->spincnts[i][j];
    }
    if(spin[1]!=0&&spin[2]!=0&&spin[3]!=0)
    {
      pbbcr->rellum[i][0]=(spin[0]+spin[2])/(spin[1]+spin[3]);
      pbbcr->rellum[i][1]=(spin[0]+spin[1])/(spin[2]+spin[3]);
      pbbcr->rellum[i][2]=(spin[0]+spin[3])/(spin[1]+spin[2]);
      pbbcr->rellum[i][3]=spin[0]/spin[3];
      pbbcr->rellum[i][4]=spin[2]/spin[3];
      pbbcr->rellum[i][5]=spin[1]/spin[3];
    }    
  }
}
void init_dat_bbc(struct dat_bbc *datp,int bx)
{
  int i;

  datp->bx=bx;

  for(i=0;i<3;i++)
  {
    datp->cnts[i]=0.;
  }
  datp->sum=0.;
}
int append_dat_bbc(struct dat_bbc *datp1,struct dat_bbc *datp2)
{
  int i;

  if(datp2->bx==datp1->bx)
  {
    for(i=0;i<3;i++)
    {
      datp2->cnts[i]+=datp1->cnts[i];
    }
    datp2->sum+=datp1->sum;
  }else
  {
    printf("datbbc2 can't be appended to datbbc1\n");
    return 0;
  }
  return 1;
}
void assign_fill(struct dat_run_bbc *pbbcr,int run)
{
  int flag;

  pbbcr->run=run;

  if(run>=12099029&&run<=12099045) flag=0;

  if(run>=12099059&&run<=12100010) flag=1;

  if(run>=12100026&&run<=12100027) flag=2;

  if(run>=12101007&&run<=12101018) flag=3;     

  if(run>=12101031&&run<=12101050) flag=4;     

  if(run>=12101060&&run<=12101064) flag=5;

  if(run>=12102030&&run<=12102054) flag=6;

  if(run>=12103009&&run<=12103029) flag=7;

  if(run>=12104003&&run<=12104019) flag=8;

  if(run>=12104085&&run<=12105016) flag=9;

  if(run>=12106017&&run<=12106055) flag=10;

  if(run>=12106067&&run<=12107014) flag=11;

  if(run>=12107026&&run<=12107048) flag=12;

  if(run>=12108004&&run<=12108020) flag=13;

  pbbcr->fill=flag;
}

void null_bx(struct dat_bbc *datp,int bx)
{
  int i;

  if(datp->bx==bx)
  {
    for(i=0;i<3;i++)
    {
      datp->cnts[i]=0;
    }
  }
}
void check_fill(struct dat_bbc *dap,int fill)
{
  int i;

  if(fill>=0&&fill<=13)
  {
    for(i=0;i<4;i++)
    {
      null_bx(dap,i);
      if(i>=0&&i<=2)
        null_bx(dap,40+i);
    }
  }
  /*
     switch(fill)
     {
     case 0:
     for(i=0;i<5;i++)
     {
     null_bx(dap,43+i);
     }
     null_bx(dap,78);
     null_bx(dap,79);
     null_bx(dap,82);
     break;
     case 1:
     null_bx(dap,78);
     null_bx(dap,79);
     break;
     case 2:
     for(i=0;i<4;i++)
     {
     null_bx(dap,i+4);
     }
     null_bx(dap,17);
     null_bx(dap,48);
     null_bx(dap,49);
     null_bx(dap,79);
     null_bx(dap,86);
     null_bx(dap,89);
     null_bx(dap,95);
     null_bx(dap,102);
     break;
     case 3:
     null_bx(dap,4);
     null_bx(dap,5);
     null_bx(dap,14);
     null_bx(dap,43);
     null_bx(dap,44);
     null_bx(dap,66);
     null_bx(dap,94);
     null_bx(dap,98);
     break;
     case 4:
     null_bx(dap,4);
     null_bx(dap,5);
     break;
     case 5:
     null_bx(dap,4);
     null_bx(dap,5);
     null_bx(dap,9);
     null_bx(dap,10);
     null_bx(dap,44);
     null_bx(dap,45);
     null_bx(dap,47);
     null_bx(dap,48);
     null_bx(dap,50);
     null_bx(dap,64);
     null_bx(dap,65);
     null_bx(dap,68);
     null_bx(dap,81);
     null_bx(dap,82);
     null_bx(dap,83);
     null_bx(dap,85);
     null_bx(dap,89);
     null_bx(dap,92);
     null_bx(dap,93);
     null_bx(dap,95);
     null_bx(dap,101);
     break;
     case 6:
     null_bx(dap,20);
     null_bx(dap,52);
     null_bx(dap,53);
  null_bx(dap,55);
  null_bx(dap,58);
  null_bx(dap,71);
  null_bx(dap,76);
  null_bx(dap,79);
  null_bx(dap,95);
  null_bx(dap,103);
  null_bx(dap,106);
  break;
  case 7:
  //     null_bx(dap,78);
  //null_bx(dap,79);
  break;
  case 8:
  break;
  case 9:
  null_bx(dap,18);
  null_bx(dap,56);
  null_bx(dap,61);
  null_bx(dap,65);
  null_bx(dap,66);
  null_bx(dap,72);
  break;
  case 10:
  null_bx(dap,10);
  null_bx(dap,11);
  null_bx(dap,12);
  null_bx(dap,13);
  null_bx(dap,21);
  null_bx(dap,25);
  null_bx(dap,43);
  null_bx(dap,44);
  null_bx(dap,78);
  null_bx(dap,79);
  null_bx(dap,90);
  null_bx(dap,97);
  null_bx(dap,98);
  break;
  case 11:
  null_bx(dap,16);
  null_bx(dap,17);
  null_bx(dap,19);
  null_bx(dap,20);
  null_bx(dap,43);
  null_bx(dap,44);
  null_bx(dap,60);
  null_bx(dap,61);
  null_bx(dap,63);
  null_bx(dap,65);
  null_bx(dap,79);
  null_bx(dap,96);
  break;
  case 12:
  null_bx(dap,22);
  null_bx(dap,67);
  null_bx(dap,68);
  null_bx(dap,70);
  null_bx(dap,78);
  null_bx(dap,79);
  break;
  case 13:
  null_bx(dap,9);
  null_bx(dap,11);
  null_bx(dap,14);
  null_bx(dap,67);
  null_bx(dap,68);
  null_bx(dap,70);
  null_bx(dap,81);
  null_bx(dap,82);
  null_bx(dap,88);
  null_bx(dap,90);
  null_bx(dap,91);
  break;
}

}
*/
