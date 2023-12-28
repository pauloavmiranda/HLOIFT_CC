
# Supplementary materials of the paper "Image Segmentation by Hierarchical Layered Oriented Image Foresting Transform Subject to Closeness Constraints" submitted to the DGMM24 conference:


## Abstract

In this work, we address the problem of image segmentation, subject to high-level constraints expected for the objects of interest. More specifically, we define closeness constraints to be used in conjunction with geometric constraints of inclusion in the Hierarchical Layered Oriented Image Foresting Transform (HLOIFT) algorithm. The proposed method can handle the segmentation of a hierarchy of objects with nested boundaries, each with its own expected boundary polarity constraint, making it possible to control the maximum distances (in a geodesic sense) between the successive nested boundaries. The method is demonstrated in the segmentation of nested objects in colored images with superior accuracy compared to its precursor methods and also when compared to some recent click-based methods.

## Authors

- Luiz F.D. Santos
- Felipe A.S. Kleine
- Paulo A.V. Miranda


## Source code

Our source code is available in the folder **"HLOIFT_closeness"**.
The code was implemented in C/C++ language, compiled with gcc 9.4.0, and tested on a Linux operating system (Ubuntu 20.04.5 LTS 64-bit), running on an Intel® Core™ i5-10210U CPU @ 1.60GHz × 8 machine. 
The code natively only supports images in the PPM format (_Portable Pixel Map_).
The [ImageMagick](https://imagemagick.org/) command-line tools are required to convert images from different formats.


The ground truth databases of public images used to evaluate our method are available in the folder **"databases"**. 
We performed quantitative experiments for colored images of 640×480 pixels divided into three databases (DB1, DB2 and DB3) with 30, 20 and 16 images, respectively, representing flat objects with nested boundaries under different viewpoints.

For each image, four segmentations were performed for different locations of the internal seed of the child object, leading to a total of 264 single-click segmentations per method. All commands executed in the experiments to generate the results, including the seeds used, are available in the subfolder **"HLOIFT_closeness/exp"**.

To compile the program, enter the folder **"HLOIFT_closeness"** and type **"make"**.
To segment an image, you must run the **"HLOIFT_closeness"** executable inside the folder with the same name **"HLOIFT_closeness"**.

### usage:

```
HLOIFT_closeness <image> <L> <radius> <seed_x> <seed_y> <hierarchy> [ground truth]
image ........... the color input image in the PPM format,
L ............... the closeness parameter L
                  (negative values indicate unconstrained results),
radius .......... the minimal distance parameter (rho),
seed_x .......... the x-coordinate of the internal object seed,
seed_y .......... the y-coordinate of the internal object seed,
hierarchy ....... hierarchy configuration file, with the boundary
                  polarity parameters (alpha_i),
ground truth .... ground truth in the PGM format.
```

As output, the program generates the label image of the resulting segmentation in file **"label.pgm"** (and its colored version in file **"label.ppm"**).

### Program execution examples:

#### To execute the HLOIFT-CC method:

The following command computes the segmentation by HLOIFT-CC for the first image in database DB1.
```
./HLOIFT_closeness ../databases/DB1/img001.ppm 40 1.5 300 200 ./txt/gsc/hierarchy_2.txt ../databases/DB1/img001_gt.pgm
```

The following command computes the segmentation by HLOIFT-CC for the first image in database DB2.
```
./HLOIFT_closeness ../databases/DB2/img001.ppm 40 1.5 465 160 ./txt/gsc/hierarchy.txt ../databases/DB2/img001_gt.pgm
```

The following command computes the segmentation by HLOIFT-CC for the first image in database DB3.
```
./HLOIFT_closeness ../databases/DB3/img001.ppm 40 1.5 234 398 ./txt/gsc/hierarchy_2.txt ../databases/DB3/img001_gt.pgm
```

#### To execute the HLOIFT-GSC method:

The following command computes the segmentation by HLOIFT-GSC for the first image in database DB1.
```
./HLOIFT_closeness ../databases/DB1/img001.ppm -1 1.5 300 200 ./txt/gsc/hierarchy_2.txt ../databases/DB1/img001_gt.pgm
```

The following command computes the segmentation by HLOIFT-GSC for the first image in database DB2.
```
./HLOIFT_closeness ../databases/DB2/img001.ppm -1 1.5 465 160 ./txt/gsc/hierarchy.txt ../databases/DB2/img001_gt.pgm
```

The following command computes the segmentation by HLOIFT-GSC for the first image in database DB3.
```
./HLOIFT_closeness ../databases/DB3/img001.ppm -1 1.5 234 398 ./txt/gsc/hierarchy_2.txt ../databases/DB3/img001_gt.pgm
```

#### To execute the regular HLOIFT method:

The following command computes the segmentation by HLOIFT for the first image in database DB1.
```
./HLOIFT_closeness ../databases/DB1/img001.ppm -1 1.5 300 200 ./txt/hierarchy_2.txt ../databases/DB1/img001_gt.pgm
```

The following command computes the segmentation by HLOIFT for the first image in database DB2.
```
./HLOIFT_closeness ../databases/DB2/img001.ppm -1 1.5 465 160 ./txt/hierarchy.txt ../databases/DB2/img001_gt.pgm
```

The following command computes the segmentation by HLOIFT for the first image in database DB3.
```
./HLOIFT_closeness ../databases/DB3/img001.ppm -1 1.5 234 398 ./txt/hierarchy_2.txt ../databases/DB3/img001_gt.pgm
```

All the above commands were taken from the files available in the subfolder **"HLOIFT_closeness/exp"**.
To run the complete experiment for a given method and database, simply run the program **"runexp"**.
The experiment output report, as a spreadsheet, will be produced in CSV format in file **"HLOIFT_closeness/out/report.csv"**.
For example, to run the experiment for method HLOIFT-CC on database DB1, simply execute:

```
./runexp exp/DB1/exp_hloiftcc.txt
```


## Contact

If you have any doubts, questions or suggestions to improve this code, please contact me at:
**pmiranda@ime.usp.br**





