
//#define APPDEBUG  1

#include "gft.h"

#include "image_io.h"
#include "hloift.h"


int *CreateSeeds(gft::sImage32 *label, int x, int y){
  int *S = NULL;
  int i,j,k,p,n;

  gft::Image32::Set(label, NIL);
  n = 1 + 2*label->ncols + 2*label->nrows - 4 + 1;
  S = (int *)calloc(n, sizeof(int));
  if(S == NULL){
    printf("Insufficient memory\n");
    exit(1);
  }
  S[0] = n;
  S[1] = x + y*label->ncols;
  label->array[y][x] = 2;
  k = 2;
  for(j = 0; j < label->ncols; j++){
    p = j;
    S[k] = p;
    label->data[p] = 0;
    k++;

    p = j + (label->nrows-1)*label->ncols;
    S[k] = p;
    label->data[p] = 0;
    k++;
  }
  for(i = 1; i < label->nrows-1; i++){
    p = i*label->ncols;
    S[k] = p;
    label->data[p] = 0;
    k++;

    p = i*label->ncols + (label->ncols-1);
    S[k] = p;
    label->data[p] = 0;
    k++;
  }

  S[0] = k-1;
  return S;
}



void AppendReportInfo(double totaltime,
		      gft::sImage32 *gtruth_label,
		      gft::sImage32 *label,
		      char *file_report){
  gft::sImage32 *bin, *gtruth_bin;
  float dice;
  int l,Lmax;
  FILE *fp;
  fp = fopen(file_report, "a");
  if(fp == NULL) return;
  Lmax = gft::Image32::GetMaxVal(gtruth_label);
  for(l = Lmax; l > 0 ; l--){
    bin        = gft::Image32::Threshold(label,        l, INT_MAX);
    gtruth_bin = gft::Image32::Threshold(gtruth_label, l, INT_MAX);
    dice = gft::Image32::DiceSimilarity(bin, gtruth_bin);
    printf("Obj: %d => Dice: %f\n", Lmax-l+1, dice);
    fprintf(fp," %f; ",dice);
    gft::Image32::Destroy(&bin);
    gft::Image32::Destroy(&gtruth_bin);
  }
  fprintf(fp," %lf\n",totaltime);
  fclose(fp);
}



int main(int argc, char **argv){
  char hierarchy[512];
  gft::sImageGraph *ig;
  gft::sLayeredGraph *lg;
  gft::sCImage *cimg = NULL, *clabel = NULL;
  gft::sCImage *lab = NULL;
  gft::sImage32 *label = NULL, *tmp = NULL, *P = NULL, *gtruth = NULL;
  float radius;
  int *S;
  int x,y,L,Wmax;
  int *hr;
  //clock_t start, end;
  struct timeval tic,toc;
  double totaltime;
  //int Sobj[2];
  
  // check number of parameters
  if(argc < 7){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"HLOIFT_closeness <image> <L> <radius> <seed_x> <seed_y> <hierarchy> [ground truth]\n");
    fprintf(stdout,"image ........... the color input image in the PPM format,\n");
    fprintf(stdout,"L ............... the closeness parameter L\n");
    fprintf(stdout,"                  (negative values indicate unconstrained results),\n");
    fprintf(stdout,"radius .......... the minimal distance parameter (rho),\n");
    fprintf(stdout,"seed_x .......... the x-coordinate of the internal object seed,\n");
    fprintf(stdout,"seed_y .......... the y-coordinate of the internal object seed,\n");
    fprintf(stdout,"hierarchy ....... hierarchy configuration file, with the boundary\n");
    fprintf(stdout,"                  polarity parameters (alpha_i),\n");
    fprintf(stdout,"ground truth .... ground truth in the PGM format.\n");
    exit(0);
  }
  
  cimg = ReadAnyCImage(argv[1]);

  L = atoi(argv[2]);
  if(L>=0 && L%2 == 0) L++;

  radius = atof(argv[3]);
  x = atoi(argv[4]);
  y = atoi(argv[5]);
  strcpy(hierarchy, argv[6]);

  if(argc >= 8){
    gtruth = ReadAnyImage(argv[7]);
  }
  
  //----------------------------
  lab = gft::CImage::RGB2Lab(cimg);
  ig = gft::ImageGraph::ByEuclideanDistance(lab, 1.5);
  label = gft::Image32::Create(cimg->C[0]);
  S = CreateSeeds(label, x, y);

#ifdef APPDEBUG
  tmp = gft::Image32::Add(label, 1);
  gft::Image32::Write(tmp, (char *)"seeds.pgm");
  gft::Image32::Destroy(&tmp);
#endif

  //------------------------------------------------
  gettimeofday(&tic,NULL);
  //start = clock();
  //------------------------------------------------

  //Sobj[0] = 1;
  //Sobj[1] = x + y*ig->ncols;
  //P = gft::ift::SC_Pred_fsum(ig, Sobj, 0.1);
  
  lg = HL_OIFT_CreateGraph(ig,
			   cimg,
			   NULL,
			   radius,
			   (char *)hierarchy,
			   &hr,
			   &Wmax,
			   &P,
			   x, y);

  if(L < 0)
    HL_OIFT_Segmentation(lg, Wmax, hr, S, label);
  else if(P != NULL)
    HL_OIFT_Closeness_Segmentation(lg, Wmax, hr, S, L, P, label, ig->A);
  else{
    printf("Invalid parameter setting\n");
    exit(1);
  }

  //------------------------------------------------  
  //end = clock();
  gettimeofday(&toc,NULL);

  //totaltime = 1000.0*(((double)(end - start))/CLOCKS_PER_SEC);
  totaltime = ((toc.tv_sec-tic.tv_sec)*1000.0 + 
	       (toc.tv_usec-tic.tv_usec)*0.001);
  printf("Total time: %f ms\n",totaltime);
  //------------------------------------------------
  
  gft::Image32::Write(label, (char *)"label.pgm");

  clabel = gft::CImage::RandomColorize(label);
  gft::CImage::Write(clabel, (char *)"label.ppm");
  gft::CImage::Destroy(&clabel);

  if(gtruth != NULL)
    AppendReportInfo(totaltime, gtruth, label, (char*)"./out/report.csv");

  //----------------------------
  gft::LayeredGraph::Destroy(&lg);
  gft::FreeIntArray(&hr);
  free(S);
  if(P != NULL)
    gft::Image32::Destroy(&P);
  gft::Image32::Destroy(&label);
  gft::ImageGraph::Destroy(&ig);
  gft::CImage::Destroy(&lab);
  gft::CImage::Destroy(&cimg);
  if(gtruth != NULL)
    gft::Image32::Destroy(&gtruth);
    
  return 0;
}

