#include<stdio.h>
#include<stdlib.h>

#include"./myscaler.h"
/*
   struct sca_val{
   int bx;
   unsigned long long valBBC[8];
   unsigned long long valZDC[8];
   unsigned long long valVPD[4];
   unsigned long long valSum;
   };
   struct sca_cnts{
   int bx;
   int spin;
   double cnts[3];

//  double cntsZDC[3];
//  double cntsZDC[3];

double sum;
};
*/
/*struct cnts_run{
  double spincnts[i][j];
  }*/
int main(int argc, char *argv[])
{
  if(argc!=3)
  {
    printf("the command line is *.exe scaler.dat spinbit.dat\n");
    return 1;
  }

  FILE *infpsca,*infpspin;
  FILE *outfbbc[3],*outfzdc[3],*outfvpd[3];


  struct sca_val scaler; // { bx, valBBC[8], valZDC[8], valVPD[4], valSum }

  // sca_cnts arrays: 0=counts 1=accidentals 2=multiples
  struct sca_cnts cntsbbc[3]; // { bx, spin cnts[3], sum }
  struct sca_cnts cntszdc[3]; // { bx, spin cnts[3], sum }
  struct sca_cnts cntsvpd[3]; // { bx, spin cnts[3], sum }
                                       // cnts[0] = east
                                       // cnts[1] = west
                                       // cnts[2] = coincidence

  char filebbc[3][100];
  char filezdc[3][100];
  char filevpd[3][100];

  char *corr[3]={"orig","accd","mult"};

  int spinbit[120];
  int i;

  // open dat file (with 22 columns)
  if((infpsca=fopen(argv[1],"r"))==NULL)
  {
    printf("can't open %s\n",argv[1]);
    return 1;
  }

  // open spinbit file
  if((infpspin=fopen(argv[2],"r"))==NULL)
  {
    printf("can't open %s\n",argv[2]);
    return 1;
  }

  // set up output files for writing
  for(i=0;i<3;i++)
  {
    sprintf(filebbc[i],"bbc_%s_%s",corr[i],argv[1]);
    sprintf(filezdc[i],"zdc_%s_%s",corr[i],argv[1]);
    sprintf(filevpd[i],"vpd_%s_%s",corr[i],argv[1]);

    outfbbc[i]=fopen(filebbc[i],"w");
    outfzdc[i]=fopen(filezdc[i],"w");
    outfvpd[i]=fopen(filevpd[i],"w");
  }

  assign_spin(infpspin,spinbit);

  fclose(infpspin);

  //  printf("%d\n",spinbit[5]);

  // read each line of data file into "scaler", which is a struct
  // containing bx, valBBC[8], valZDC[8], valVPD[4], valSum
  while(scan_sca_val(infpsca,&scaler)==22)
  {
    // set bx and spinbit for sca_cnts structs (bbc, zdc, vpd)
    // -- note that there are three of each of these structs
    for(i=0;i<3;i++)
    {
      cntsbbc[i].bx=scaler.bx;
      //      printf("%d\n",scaler.bx);
      cntsbbc[i].spin=spinbit[scaler.bx];
      //      printf("%d\n",cntsbbc.spin);
      cntszdc[i].bx=scaler.bx;
      cntszdc[i].spin=spinbit[scaler.bx];
      cntsvpd[i].bx=scaler.bx;
      cntsvpd[i].spin=spinbit[scaler.bx];
    }

    cnts_bbc(&scaler,&cntsbbc[0]);
    accidental(&cntsbbc[1],&cntsbbc[0]);
    multiple(&cntsbbc[2],&cntsbbc[1]);

    cnts_zdc(&scaler,&cntszdc[0]);
    accidental(&cntszdc[1],&cntszdc[0]);
    multiple(&cntszdc[2],&cntszdc[1]);

    cnts_vpd(&scaler,&cntsvpd[0]);
    accidental(&cntsvpd[1],&cntsvpd[0]);
    multiple(&cntsvpd[2],&cntsvpd[1]);


    for(i=0;i<3;i++)
    {
      print_sca_cnts(outfbbc[i],&cntsbbc[i]);	// prints columns [bx][spin][cnts][acc][mult][sum]  
      print_sca_cnts(outfzdc[i],&cntszdc[i]);
      print_sca_cnts(outfvpd[i],&cntsvpd[i]);
    }
  }

  fclose(infpsca);

  for(i=0;i<3;i++)
  {
    fclose(outfzdc[i]);
  }

  return 0;
}



