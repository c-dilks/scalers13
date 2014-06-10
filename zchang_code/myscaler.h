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
  /*                                                                            
  double cntsZDC[3];                                                            
  double cntsZDC[3];                                                            
  */
  double sum;
};
struct cnts_rel{
  int run;

  double spincnts[3][4];
  double sum;
  double rellum[3][6];
};
