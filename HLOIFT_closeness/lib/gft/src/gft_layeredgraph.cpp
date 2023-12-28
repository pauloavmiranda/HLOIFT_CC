
#include "gft_layeredgraph.h"

namespace gft{
  namespace LayeredGraph{

    sLayeredGraph *Create(int nlayers, int nnodesperlayer){
      sLayeredGraph *lg;
      lg = (sLayeredGraph *)calloc(1, sizeof(sLayeredGraph));
      if(lg == NULL)
	gft::Error((char *)MSG1,(char *)"LayeredGraph::Create");
      lg->nnodesperlayer = nnodesperlayer;
      lg->nlayers = nlayers;
      lg->graph = gft::Graph::Create(nlayers*nnodesperlayer, 8+8*nlayers, NULL);
      return lg;
    }

    
    void Destroy(sLayeredGraph **lg){
      sLayeredGraph *tmp;
      if(lg == NULL) return;
      tmp = *lg;
      if(tmp != NULL){
	gft::Graph::Destroy(&tmp->graph);
	free(tmp);
      }
      *lg = NULL;
    }

    //2D images:
    void SetArcs(sLayeredGraph *lg, sImageGraph *sg, int l){
      int p, q, i, n, src, dest;
      int p_x, p_y, q_x, q_y;
      float wpq;
      n = sg->ncols*sg->nrows;
      if(n != lg->nnodesperlayer)
	gft::Error((char *)"Incompatible number of nodes per layer.",
		   (char *)"SetArcs");
      for(p = 0; p < n; p++){
	p_x = p % sg->ncols;
	p_y = p / sg->ncols;
	for(i = 1; i < sg->A->n; i++){
	  q_x = p_x + sg->A->dx[i];
	  q_y = p_y + sg->A->dy[i];

          if ((q_x >= 0)&&(q_x < sg->ncols)&&
	      (q_y >= 0)&&(q_y < sg->nrows)){
	    
	    q = q_x + q_y*sg->ncols;
	    
	    wpq = (float)((sg->n_link[p])[i]);
	    //printf("Wpq %f\n",wpq);
	    
	    src = p + l*n;
	    dest = q + l*n;
	    gft::Graph::AddDirectedEdge(lg->graph, src, dest, wpq);
	  }
	  else{ 
	    src = p + l*n;
	    gft::Graph::AddDirectedEdge(lg->graph, src, src, NIL);
	  }
	}
      }
    }


    void SetArcs(sLayeredGraph *lg,
		 int l_orig, int l_dest,
		 int ncols,
		 float w, float r){
      SetArcs(lg, l_orig, l_dest,
	      ncols, w, r, 1.5);
    }

    
    void SetArcs(sLayeredGraph *lg,
		 int l_orig, int l_dest,
		 int ncols,
		 float w, float r, float r_neighborhood){
      sAdjRel *A;
      int n, p, q, i, src, dest;
      int p_x, p_y, q_x, q_y;
      int nrows;
      nrows = lg->nnodesperlayer/ncols;
      n = lg->nnodesperlayer;
      A = gft::AdjRel::Circular(r);
      //printf("A-n old: %d\n",A->n);
      gft::AdjRel::Reduce2Boundary(&A, r_neighborhood);
      //printf("A-n new: %d\n",A->n);
      
      //printf("Entro set Arcs 2 \n");

      for(p = 0; p < n; p++){
	p_x = p % ncols;
	p_y = p / ncols;

	for(i = 0; i < A->n; i++){ //i = 1;
	  q_x = p_x + A->dx[i];
	  q_y = p_y + A->dy[i];
          //printf("p_x, p_y = %d, %d \n",p_x,p_y);
          //printf("q_x, q_y = %d, %d \n",q_x,q_y);

	  if ((q_x >= 0)&&(q_x < ncols)&&
	      (q_y >= 0)&&(q_y < nrows)){
	    
	    q = q_x + q_y*ncols;
	    
	    src  = p + l_orig*n;
	    dest = q + l_dest*n;
	    
	    //printf("Wpq_L %f\n",w);
	    gft::Graph::AddDirectedEdge(lg->graph, src, dest, w);
	  }
	  /*
	  else{ 
	  q = q_x + q_y*ncols;
	  
	  src  = p + l_orig*n;
	  dest = q + l_dest*n;
	  gft::Graph::AddDirectedEdge(lg->graph, src, dest, NIL);
	  }
	  */
	}
      }
      gft::AdjRel::Destroy(&A);
    }



    //ND images:
    void SetArcs(sLayeredGraph *lg, sGraph *g, int l){
      int p,q,i,w,src,dest,n;
      n = lg->nnodesperlayer;
      for(p = 0; p < g->nnodes; p++){
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];
	  w = g->nodes[p].Warcs[i];
	  src = p + l*n;
	  dest = q + l*n;
	  gft::Graph::AddDirectedEdge(lg->graph, src, dest, w);
	}
      }
    }
    

    void SetArcs(sLayeredGraph *lg,
		 int l_orig, int l_dest,
		 float w){
      int p,q,n,i,src,dest;
      n = lg->nnodesperlayer;
      for(p = 0; p < n; p++){
	src = p + l_orig*n;
	for(i = 0; i < lg->graph->nodes[src].outdegree; i++){
	  q = lg->graph->nodes[src].adjList[i];
	  if(q >= l_orig*n && q < (l_orig+1)*n){
	    dest = q - l_orig*n + l_dest*n;
	    gft::Graph::AddDirectedEdge(lg->graph, src, dest, w);
	  }
	}
	dest = p + l_dest*n;
	gft::Graph::AddDirectedEdge(lg->graph, src, dest, w);      
      }
    }
    
    
    void TransposeLayer(sLayeredGraph *lg, int l){
      int p, q, n, i;
      float wpq, wqp;
      n = lg->nnodesperlayer;
      for(p = l*n; p < (l+1)*n; p++){
	for(i = 0; i < lg->graph->nodes[p].outdegree; i++){
	  q = lg->graph->nodes[p].adjList[i];
	  if(q >= l*n && q < (l+1)*n){
	    if(q > p){
	      wpq = lg->graph->nodes[p].Warcs[i];
	      wqp = gft::Graph::GetArcWeight(lg->graph, q, p);
	      gft::Graph::UpdateDirectedEdge(lg->graph, q, p, wpq);
	      gft::Graph::UpdateDirectedEdge(lg->graph, p, q, wqp);
	    }
	  }
	}
      }
    }



    void RemoveCuttingEdges(sLayeredGraph *lg, int l,
			    int *label, int lb){
      int p, q, n, i;
      n = lg->nnodesperlayer;
      for(p = l*n; p < (l+1)*n; p++){
	for(i = 0; i < lg->graph->nodes[p].outdegree; i++){
	  q = lg->graph->nodes[p].adjList[i];
	  if(q >= l*n && q < (l+1)*n){
	    if(label[p] == lb && label[q] != lb){
	      gft::Graph::RemoveEdge(lg->graph, p, q);
	    }
	  }
	}
      }
    }

    
    
  } //end LayeredGraph namespace
} //end gft namespace

