/* Chris Perkins - <cwperkins@lbl.gov> - 2011115
 *
 * Example reader for ScalerII Histogrammed Data Files
 *
 */

#include <stdio.h>
#include <fcntl.h>
#include <string.h>

struct hist_hdr {
  unsigned int version;
  unsigned int bdnum;
  unsigned int runnum;
  unsigned int count;
  unsigned int num_words;
  unsigned int rs_count_lo;
  unsigned int rs_count_hi;
  unsigned int status1;
  unsigned int status2;
  unsigned int status3;
};

// recursive binary number printer
void bin(unsigned n)
{
  if(n>1) bin(n/2);
  printf("%d",n%2);
};

int main(int argc, char *argv[]) {
  struct hist_hdr hist_hdr_l;
  char filename_l[100];
  int fd;
  int bytes_read;
  int i;
  int ki;
  unsigned int chn_cnt[3];
  unsigned long long channel;
  unsigned long long count;
  unsigned long long sum = 0;
  int debug_l = 0; // debug flag

  printf("Usage: scaler2_reader [filename] [debug]\n");

  if(argc > 1) strcpy(filename_l, argv[1]);
  else sprintf(filename_l, "out.dat");
  if(argc > 2) debug_l = atoi(argv[2]);
  printf("Opening datafile: %s\n", filename_l);

  // Open File
  fd = open(filename_l, O_RDONLY);
  if(fd == -1) {
    printf("Error opening datafile: %s\n", filename_l);
    return -1;
  }

  // Read Version
  bytes_read = read(fd, &hist_hdr_l, 4);
  // Data Format Version 1
  if(hist_hdr_l.version == 1) 
  {  
    // Read Rest of Header
    bytes_read = read(fd, &(hist_hdr_l.bdnum), sizeof(struct hist_hdr) - 4);
    if(bytes_read != (sizeof(struct hist_hdr) - 4)) {
      printf("Error reading header\n");
      close(fd);
      return -1;
    }

    // print hist_hdr_l (read buffer --> struct)
    if(debug_l==2)
    {
      printf("\n");
      printf("[+] BEGIN HEADER BUFFER [+]\n");
      printf(" version=%d\n",hist_hdr_l.version);
      printf(" bdnum=%d\n",hist_hdr_l.bdnum);
      printf(" runnum=%d\n",hist_hdr_l.runnum);
      printf(" count=%d\n",hist_hdr_l.count);
      printf(" num_words=%d\n",hist_hdr_l.num_words);
      printf(" rs_count_lo=%d\n",hist_hdr_l.rs_count_lo);
      printf(" rs_count_hi=%d\n",hist_hdr_l.rs_count_hi);
      printf(" status1=%d\n",hist_hdr_l.status1);
      printf(" status2=%d\n",hist_hdr_l.status2);
      printf(" status3=%d\n",hist_hdr_l.status3);
      printf("[+] END HEADER BUFFER [+]\n");
      printf("\n");
    };

    // Print Header
    printf("Runnumber = %8.8d, BoardNumber = %2.2d\n", hist_hdr_l.runnum, hist_hdr_l.bdnum);
    printf("  Data Format Version = %d, Count = %d\n", hist_hdr_l.version, hist_hdr_l.count);
    printf("  Num Data Words = %d (0x%8.8x)\n", hist_hdr_l.num_words, hist_hdr_l.num_words);
    printf("  Num Channels = %d (0x%8.8x),  RS Count = 0x%8.8x%8.8x\n", hist_hdr_l.num_words / 3, hist_hdr_l.num_words / 3, hist_hdr_l.rs_count_hi, hist_hdr_l.rs_count_lo);
    printf("  Status Words = 0x%8.8x  0x%8.8x  0x%8.8x\n", hist_hdr_l.status1, hist_hdr_l.status2, hist_hdr_l.status3);

    // Read Channels and Counts
    for(i=0; i < (hist_hdr_l.num_words / 3) ; i++) {
      bytes_read = read(fd, chn_cnt, 12);
      if(bytes_read != 12) {
        printf("Error reading Channel/Count\n");
        close(fd);
        return -1;
      }
      channel = ((unsigned long long)(chn_cnt[2]) << 16) | ((unsigned long long)(chn_cnt[1]) >> 16);
      count = (((unsigned long long)(chn_cnt[1]) & 0xffff) << 32) | (unsigned long long)(chn_cnt[0]);
      sum += count;
      if(debug_l==1) printf("Channel 0x%12.12llx = 0x%12.12llx Counts\n", channel, count);
      if(debug_l==3) 
      {
        printf("chn_cnt: ");
        for(ki=0; ki<3; ki++) { bin(chn_cnt[ki]); printf("  "); };
        printf("\n");
      };
    }
    printf("\n");
    printf("RS Count = 0x%8.8x%8.8x, Count Sum = 0x%16.16llx\n", hist_hdr_l.rs_count_hi, hist_hdr_l.rs_count_lo, sum);
  } 
  else
  {
    printf("Unknown Data File Format Version: %d\n", hist_hdr_l.version);
  }

  close(fd);

  return 0;
}

