
#include "image_io.h"


gft::sCImage *ReadAnyCImage(char *file){
  gft::sCImage *cimg;
  char command[512];
  int s;

  s = strlen(file);
  if(strcasecmp(&file[s-3], "ppm") == 0){
    cimg = gft::CImage::Read(file);
  }
  else{
    sprintf(command, "convert %s cimage_tmp.ppm", file);
    system(command);
    cimg = gft::CImage::Read("cimage_tmp.ppm");
    system("rm cimage_tmp.ppm");
  }
  return cimg;
}


gft::sImage32 *ReadAnyImage(char *file){
  gft::sImage32 *img;
  char command[512];
  int s;

  s = strlen(file);
  if(strcasecmp(&file[s-3], "pgm") == 0){
    img = gft::Image32::Read(file);
  }
  else{
    sprintf(command, "convert %s image_tmp.pgm", file);
    system(command);
    img = gft::Image32::Read("image_tmp.pgm");
    system("rm image_tmp.pgm");
  }

  return img;
}


