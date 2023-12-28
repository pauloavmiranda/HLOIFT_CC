
#ifndef _GFT_LAYEREDGRAPH_H_
#define _GFT_LAYEREDGRAPH_H_

#include "gft_common.h"
#include "gft_stack.h"
#include "gft_graph.h"
#include "gft_imagegraph.h"


namespace gft{
  namespace LayeredGraph{

    struct sLayeredGraph {
      int nlayers;
      int nnodesperlayer;
      sGraph *graph;
    };

    sLayeredGraph *Create(int nlayers, int nnodesperlayer);
    void Destroy(sLayeredGraph **lg);

    //2D images:
    void SetArcs(sLayeredGraph *lg, sImageGraph *sg, int l);
    void SetArcs(sLayeredGraph *lg,
		 int l_orig, int l_dest,
		 int ncols,
		 float w, float r);
    void SetArcs(sLayeredGraph *lg,
		 int l_orig, int l_dest,
		 int ncols,
		 float w, float r, float r_neighborhood);

    //ND images:
    void SetArcs(sLayeredGraph *lg, sGraph *g, int l);

    void SetArcs(sLayeredGraph *lg,
		 int l_orig, int l_dest,
		 float w);
    
    void TransposeLayer(sLayeredGraph *lg, int l);

    void RemoveCuttingEdges(sLayeredGraph *lg, int l,
			    int *label, int lb);

    
  } //end LayeredGraph namespace

  typedef LayeredGraph::sLayeredGraph sLayeredGraph;

} //end gft namespace

    
#endif


    
