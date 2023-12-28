
#include "gft.h"


int main(int argc, char **argv){
  FILE *fp;
  //struct timeval tic,toc;
  //double totaltime;
  char command[512];
  int ret, ntests;
  
  // check number of parameters
  if(argc < 2){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"runexp <exp_file>\n");
    exit(0);
  }

  fp = fopen(argv[1], "r");
  if(fp == NULL){
    printf("File not found.\n");
    exit(1);
  }

  system("rm ./out/report.csv -f");
  system("touch ./out/report.csv");
  system("echo \"Dice obj1; Dice obj2; Time\" > ./out/report.csv");
  
  //gettimeofday(&tic,NULL);
  //------------------------------------------------
  ntests = 0;
  while(1){
    ret = fscanf(fp," %[^\n]", command);
    if(ret==EOF)
      break;
    gft::String::Trim(command);

    //Ignore a comment line:
    if(command[0]=='#')
      continue;    

    ntests++;
    printf("***%s***\n",command);
    system(command);
  }
  
  //------------------------------------------------  
  //gettimeofday(&toc,NULL);

  //totaltime = ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);
  //printf("Total time: %f ms\n",totaltime);

  printf("Number of tests: %d\n", ntests);
  
  fclose(fp);
  
  return 0;
}
