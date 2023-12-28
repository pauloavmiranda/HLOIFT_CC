
#ifndef _GFT_GRAPH_H_
#define _GFT_GRAPH_H_

#include "gft_common.h"
#include "gft_stack.h"
#include "gft_image32.h"
#include "gft_filtering.h"
#include "gft_imagegraph.h"
#include "gft_bmap.h"
#include "gft_cimage32f.h"
#include "gft_adjrel.h"
#include "gft_draw.h"
#include "gft_analysis.h"


namespace gft{

  /**
   * \brief Common definitions and functions to manipulate a generic graph.
   */
  namespace Graph{

    /**
     * \brief Structure representing a node of a generic graph.
     */
    struct sGraphNode {
      int outdegree; /**< The current number of adjacent nodes for the given node.*/
      int arraysize; /**< The actual total allocated length of the array used to store the adjacency list.*/
      int id;
      int *adjList; /**< Adjacency list (list of nodes that are adjacent to the given node) implemented as a dynamically allocated array with current length in use given by outdegree.*/
      int *Warcs; /**< The corresponding arc weights for the adjacency list.*/
    };


    /**
     * \brief Structure representing a generic graph.
     */
    struct sGraph {
      int type;
      int nnodes; /**< The total number of nodes in the graph.*/
      sGraphNode *nodes; /**< The array of nodes in the graph, with indexes in the range [0,  nnodes-1].*/
      int        *Wnodes;
    };


    /**
     * \brief A constructor.
     *
     * A utility function that creates a graph of 'nnodes' vertices.
     */
    sGraph *Create(int nnodes, int outdegree, int *Wnodes);

    /**
     * \brief A copy constructor.
     */
    sGraph *Clone(sGraph *graph);
    /**
     * \brief A copy constructor.
     */
    sGraph *Clone(sImageGraph *sg);

    /**
     * \brief Updates the weights of the arcs to incorporate the boundary polarity constraint.
     *
     * It is usually used before computing the OIFT 
     * (Oriented Image Foresting Transform).
     */
    void Orient2Digraph(sGraph *graph, 
			int *data,
			float per);

    /**
     * \brief Updates the weights of the arcs to incorporate the Geodesic Star Convexity constraint.
     *
     * It is usually used before computing the OIFT 
     * (Oriented Image Foresting Transform).
     */
    void Orient2DigraphOuter(sGraph *graph,
			     int *P_sum);
    
    /**
     * \brief Computes the transpose of a directed graph.
     *
     * Transpose of a directed graph G is another directed graph on the 
     * same set of vertices with all of the arcs reversed compared to 
     * the orientation of the corresponding arcs in G. That is, if G 
     * contains an arc (u, v) then the transpose of G contains an 
     * arc (v, u) and vice versa.
     */
    sGraph *Transpose(sGraph *graph);

    /**
     * \brief A destructor.
     */
    void   Destroy(sGraph **graph);

    /**
     * \brief Returns the total number of arcs in the graph.
     */
    int GetNumberOfArcs(sGraph *graph);

    /**
     * \brief Returns the weight of the arc (src, dest) in the graph.
     */
    int GetArcWeight(sGraph *graph, int src, int dest);
    /**
     * \brief Returns the index of a given arc (src, dest) in the adjacency list of src.
     *
     * Warning: This index may be affected by the removal of other arcs, becoming inconsistent in future operations.
     */
    int GetArcIndex(sGraph *graph, int src, int dest);
    
    /** 
     * \brief Adds an edge to an undirected graph.
     *
     * Since undirected graphs are internally represented as 
     * symmetric graphs, this function adds two arcs (src, dest) 
     * and (dest, src) to the graph. 
     */
    void AddEdge(sGraph *graph, int src, int dest, int w);
    
    /** 
     * \brief Adds an arc to a directed graph.
     */
    void AddDirectedEdge(sGraph *graph, int src, int dest, int w);
    /** 
     * \brief Adds an arc to a directed graph only if it is an unrepeated arc.
     */
    void AddDirectedEdgeIfNonexistent(sGraph *graph, int src, int dest, int w);

    /**
     * \brief Updates the weight of an edge of the undirected graph. The edge is created if it does not exist.
     */
    void UpdateEdge(sGraph *graph, int src, int dest, int w);

    /**
     * \brief Updates the weight of an arc of the graph. The arc is created if it does not exist.
     */    
    void UpdateDirectedEdge(sGraph *graph, int src, int dest, int w);

    /**
     * \brief Updates the weight of an arc of the graph only if the new value is higher. The arc is created if it does not exist.
     */
    void UpdateDirectedEdgeIfHigher(sGraph *graph, int src, int dest, int w);

    /**
     * \brief Removes an edge from the undirected graph.
     */
    void RemoveEdge(sGraph *graph, int src, int dest);

    /**
     * \brief Removes an arc from the graph.
     */
    void RemoveDirectedEdge(sGraph *graph, int src, int dest);

    /**
     * \brief Returns the value of the highest arc weight of the graph.
     */
    int GetMaximumArc(sGraph *graph);

    /**
     * \brief Returns the value of the lowest arc weight of the graph.
     */
    int GetMinimumArc(sGraph *graph);
    
    /** 
     * \brief Computes the Tarjan's algorithm.
     * 
     * Tarjan's algorithm is a procedure for finding strongly connected 
     * components of a directed graph.
     */
    int *Tarjan(sGraph *graph);
    /** 
     * \brief Computes the Tarjan's algorithm.
     * 
     * Tarjan's algorithm is a procedure for finding strongly connected 
     * components of a directed graph.
     */    
    int *Tarjan(sGraph *graph, int p);
    /** 
     * \brief Computes the Tarjan's algorithm.
     * 
     * Tarjan's algorithm is a procedure for finding strongly connected 
     * components of a directed graph.
     */
    int *Tarjan(sGraph *graph, int *V, int n);

    void ChangeType(sGraph *graph, int type);

    /**
     * \brief Replaces a given arc weight with a new value.
     */    
    void ChangeWeight(sGraph *graph,
		      int old_val, int new_val);

    /**
     * \brief Applies an increasing transformation to the arc weights to enhance their differences.
     */
    void Pow(sGraph *graph, int power, int max);

    /**
     * \brief Adds/updates the arc weights to incorporate the Hedgehog shape constraint.
     */
    void HedgehogDigraph(sGraph *graph,
			 sImage32 *cost,
			 float theta,
			 float r);

    /**
     * \brief Computes the transitive reduction of a DAG representing the Local Band constraint.
     */
    void LocalBandReduction(sGraph **graph,
			    sImage32 *cost,
			    int delta);

    /**
     * \brief Computes the transitive reduction of a DAG representing the Local Band constraint.
     */
    void LocalBandReduction(sGraph **graph,
			    int *cost,
			    int delta);

    /**
     * \brief Adds/updates the arc weights to incorporate the Local Band shape constraint.
     */
    void LocalBandConstraint(sGraph *graph,
			     sImage32 *cost,
			     int delta,
			     float r);
    
    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG(sCImage *cimg,
		sImage32 *label);
    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG_lab(sCImage *cimg,
		    sImage32 *label);

    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG_lab(sCImage *cimg,
		    sImage32 *label,
		    int c1,
		    int c2,
		    float polarity);

    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG_lab(sCImage *cimg,
		    sImage32 *label,
		    int c1,
		    int c2,
		    float polarity,
		    float R);
    
    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG_lab(sCImage *cimg,
		    sImage32 *label,
		    float polarity);

    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG_lab(sCImage *cimg,
		    sImage32 *label,
		    float polarity,
		    float R);
    
    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG(sImage32 *img,
		sImage32 *label);
    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG(sImage32 *img,
		sImage32 *label,
		float polarity);

    /**
     * \brief Computes a Region Adjacency Graph (RAG).
     */
    sGraph *RAG(sImage32 *img,
		sImage32 *label,
		float polarity,
		float R);

    /**
     * \brief Generates a drawing for viewing a provided RAG.
     */
    sCImage *DrawRAG(sGraph *G,
		     sImage32 *label, float scale,
		     float R, float r, float w, float h);

    /**
     * \brief Splits the graph into two graphs according to a labeling provided.
     */    
    sGraph *Split(sGraph **G, int *label, int lb);



    /**
     * \brief Performs contrast enhancement by linear stretch of the arc weights.
     */
    void LinearStretch(sGraph *G,
		       int f1, int f2,
		       int g1, int g2);

    /**
     * \brief Computes the weighted arithmetic mean of the corresponding arcs of two graphs.
     */
    sGraph *WeightedMean(sGraph *G1, 
			 sGraph *G2,
			 float w2);
    
  } //end Graph namespace

  typedef Graph::sGraph sGraph;
  typedef Graph::sGraphNode sGraphNode;
  
} //end gft namespace

    
#endif

