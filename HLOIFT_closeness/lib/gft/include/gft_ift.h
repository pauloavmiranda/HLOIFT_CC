
#ifndef _GFT_IFT_H_
#define _GFT_IFT_H_

#include "gft_common.h"
#include "gft_imagegraph.h"
#include "gft_pqueue32.h"
#include "gft_queue.h"
#include "gft_heap.h"
#include "gft_heap32.h"
#include "gft_heap_lex.h"
#include "gft_heap32fi_lex.h"
#include "gft_heap32fifi_lex.h"
#include "gft_heap32fiif_lex.h"
#include "gft_heap32fif_lex.h"
#include "gft_heappair.h"
#include "gft_stack.h"
#include "gft_graph.h"
#include "gft_set.h"
#include "gft_morphology.h"
#include "gft_layeredgraph.h"
#include "gft_analysis.h"
#include "gft_image32.h"

//----------------------------------------
//  Method wrappers:

//The 'label' scene should be pre-initialized as follows:
//  label->data[p]=NIL, unlabeled voxel.
//  label->data[p]=0,   background voxel.
//  label->data[p]=1,   object#1 voxel.
//  label->data[p]=2,   object#2 voxel.
//  ...
//---------------------------------------

namespace gft{
  namespace ift{


    //---------------------------------------
    /* Star Convexity - ICIP 2013 */
    sImage32 *SC_Pred_fsum(sImageGraph *sg,
			   int *S, float power);
    void SC_IFT(sImageGraph *sg,
		int *S,
		sImage32 *label,
		sImage32 *P_sum);
    //---------------------------------------
    
    //Outer Cut:
    int GetEnergy_Min(sImageGraph *sg,
		      sImage32 *label,
		      int lb);
    int GetEnergy_Max(sImageGraph *sg,
		      sImage32 *label,
		      int lb);
    long long GetEnergy_Sum(sImageGraph *sg,
			    sImage32 *label,
			    int lb);

    int GetEnergy_Min(sGraph *graph,
		      int *label,
		      int lb);

    float GetEnergy_Mean(sGraph *graph,
			 int *label,
			 int lb);
    
    //---------------------------------------
    // OIFT:
    
    //Outer Cut:
    void OIFT(sImage32 *W,
	      sAdjRel *A,
	      sImage32 *img,
	      float per,
	      int *S,
	      sImage32 *label);

    //Outer Cut (the polarity must be embedded in the graph):
    void OIFT(sImageGraph *sg,
	      int *S,
	      sImage32 *label);

    void OIFT_MinMax(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter);

    void OIFT_MaxMin(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter);
    
    void OIFT_TZ2Bkg(sImageGraph *sg,
		     int *S,
		     sImage32 *label);
    void OIFT_TZ2Obj(sImageGraph *sg,
		     int *S,
		     sImage32 *label);
    
    void OIFT_TZ(sImageGraph *sg,
		 int *S,
		 sImage32 *label);

    bool isOIFT(sImageGraph *sg,
		int *S,
		sImage32 *Slabel,
		sImage32 *label,
		sImage32 *pred,
		sImage32 *ord,
		bool complete_check);
    
    bool isOIFT_Segmentation(sImageGraph *sg,
			     int *S,
			     sImage32 *Slabel,
			     sImage32 *label);

    bool isOIFT_Forest(sImageGraph *sg,
		       int *S,
		       sImage32 *Slabel,
		       sImage32 *pred,
		       sImage32 *ord);

    bool isForest(sImageGraph *sg,
		  sImage32 *pred);

    //Remover:
    void OIFT_guided(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *ord);
    
    //Outer Cut (the polarity must be embedded in the graph):
    void OIFT(sGraph *graph,
	      sGraph *transpose,
	      int *S,
	      int *label);

    //Inner Cut (the polarity must be embedded in the graph):
    void OIFT_in(sImageGraph *sg,
		 int *S,
		 sImage32 *label);

    void OIFT_Heap(sImageGraph *sg,
		   int *S,
		   sImage32 *label);
    void OIFT_Heap(sGraph *graph,
		   sGraph *transpose,
		   int *S,
		   int *label);
        
    //Outer Cut (the polarity must be embedded in the graph):
    //Energy based implementations.

    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label,
	       int e_max);

    void EOIFT(sGraph *graph,
	       sGraph *transpose,
	       int *S,
	       int *label,
	       int e_max);

    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max);
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label,
		      float e_max);

    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max);
    
    void EOIFT_Heap_2(sGraph *graph,
		      sGraph *transpose,
		      int *S,
		      int *label,
		      float e_max);

    //---------------------------------------

    void IFT_fmax(sGraph *graph,
		  int *S,
		  int *label,
		  int *cost);
     
    void IFT_fmax_Heap(sGraph *graph,
		       int *S,
		       int *label,
		       float *cost);

    void IFT_feuc(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred);
    
    void IFT_fw(sImageGraph *g,
		int *S,
		sImage32 *label,
		sImage32 *cost,
		sImage32 *pred);

    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost);

    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost,
		int *pred);
        
    void IFT_fw_Heap(sGraph *graph,
		     int *S,
		     int *label,
		     float *cost,
		     int *pred);
    
    //---------------------------------------
    sImage32 *Cost_fmin(sImageGraph *sg,
			int *S, int lb,
			sImage32 *label);

    //---------------------------------------

    int *GetSeedsByLabel(int *S,
			 sImage32 *label,
			 int lb);
    int *GetAllInternalSeedsByLabel(int *S,
				    sImage32 *label,
				    int lb,
				    int *hr,
				    int nlayers);
    sImageGraph *GetPolarityGraph(sImageGraph *graph,
				  sCImage *cimg,
				  sImage32 *img,
				  int pol);
    
    void HL_OIFT(sImageGraph *graph,
		 sCImage *cimg,
		 sImage32 *img,
		 float radius,
		 char *hierarchy,
		 int *S,
		 sImage32 *label);
    
    void HL_OIFT(sLayeredGraph *lg,
		 int Wmax, int *S, int *L, int *hr);

    void HL_OIFT_2(sLayeredGraph *lg,
		   int Wmax, int *S, int *L, int *hr);

    //---------------------------------------

    //Compute watershed by fpeak from markers
    void IFT_fpeak(sImage32 *grad,
		   sAdjRel *A,
		   sImage32 *label);
    
    //Compute watershed by fwv from markers
    void IFT_fwv(sImage32 *grad,
		 sAdjRel *A,
		 sImage32 *label);

    //--------------------------------------
    
  } //end ift namespace
} //end gft namespace

#endif

