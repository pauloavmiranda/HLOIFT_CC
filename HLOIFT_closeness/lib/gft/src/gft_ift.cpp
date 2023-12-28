
#include "gft_ift.h"
#include <queue>

namespace gft{
  namespace ift{

    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_conquer_path(int p,
			 sImageGraph *sg,
			 sImage32 *P_sum, 
			 sImage32 *V,
			 sPQueue32 *Q,
			 sImage32 *label);
    
    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_prune_tree(int p,
		       sImageGraph *sg,
		       sImage32 *P_sum, 
		       sImage32 *V,
		       sPQueue32 *Q,
		       sQueue *Qfifo,
		       sImage32 *label);
    
    //----------------------------------


    void IFT_feuc(sImageGraph *g,
		  int *S,
		  sImage32 *label,
		  sImage32 *cost,
		  sImage32 *pred){
      sImage32 *Dx=NULL,*Dy=NULL;
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst,dx,dy;
      Pixel u,v;
      sAdjRel *A = g->A;

      Dx = Image32::Create(cost->ncols, cost->nrows);
      Dy = Image32::Create(cost->ncols, cost->nrows);
      
      n = g->ncols*g->nrows;
      Q = PQueue32::Create(2*(label->ncols+label->nrows), n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (gft::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;

	    dx  = Dx->data[p] + abs(A->dx[i]);
	    dy  = Dy->data[p] + abs(A->dy[i]);
	    cst = dx*dx + dy*dy; //(g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      Dx->data[q] = dx;
	      Dy->data[q] = dy;
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      gft::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      Image32::Destroy(&Dx);
      Image32::Destroy(&Dy);      
    }


    
    void IFT_fw(sImageGraph *g,
		int *S,
		sImage32 *label,
		sImage32 *cost,
		sImage32 *pred){
      sPQueue32 *Q=NULL;
      int i,p,q,n,cst;
      Pixel u,v;
      sAdjRel *A = g->A;
      
      n = g->ncols*g->nrows;
      Q = PQueue32::Create(g->Wmax+2, n, cost->data);

      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) cost->data[p] = INT_MAX;
	else                    cost->data[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if (gft::Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){	  
	    q = v.x + v.y*label->ncols;
	    
	    cst = (g->n_link[p])[i];
	    
	    if(cst < cost->data[q]){	    
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      cost->data[q]  = cst;
	      pred->data[q]  = p;
	      label->data[q] = label->data[p];
	      gft::PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
    }
    
    //----------------------------------
    
    int GetEnergy_Min(sImageGraph *sg,
		      sImage32 *label,
		      int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      int energy,w;

      energy = INT_MAX;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy = MIN(energy, w);
	    }
	  }
	}
      }
      return energy;
    }


    int GetEnergy_Max(sImageGraph *sg,
		      sImage32 *label,
		      int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      int energy,w;

      energy = INT_MIN;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy = MAX(energy, w);
	    }
	  }
	}
      }
      return energy;
    }
    

    long long GetEnergy_Sum(sImageGraph *sg,
			    sImage32 *label,
			    int lb){
      sAdjRel *A;
      int u_x,u_y,v_x,v_y,p,q,n,i;
      long long energy,w;

      energy = 0;
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(label->data[p] == lb && label->data[q] != lb){
	      w = (sg->n_link[p])[i];
	      energy += w;
	    }
	  }
	}
      }
      return energy;
    }



    int GetEnergy_Min(sGraph *graph,
		      int *label,
		      int lb){
      int p,q,i;
      int energy,w;

      energy = INT_MAX;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    energy = MIN(energy, w);
	  }
	}
      }
      return energy;
    }



    float GetEnergy_Mean(sGraph *graph,
			 int *label,
			 int lb){
      int p,q,i;
      int w, n = 0;
      float sum = 0.0;
      for(p = 0; p < graph->nnodes; p++){
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];
	  if(label[p] == lb && label[q] != lb){
	    sum += w;
	    n++;
	  }
	}
      }
      return (sum/n);
    }

    
    
    
    /*Weighted Distance Transform.*/
    sImage32 *SC_Pred_fsum(sImageGraph *sg,
			   int *S,
			   float power){
      sHeap *Q=NULL;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      int u_x,u_y,v_x,v_y;
      sImage32 *pred;
      sAdjRel *A;
      float *Dpq;
      
      n    = sg->ncols*sg->nrows;
      pred = Image32::Create(sg->ncols, sg->nrows);
      cost = gft::AllocFloatArray(n);
      Q = Heap::Create(n, cost);
      A = sg->A;

      //--------------------
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
      }
      //--------------------
      
      Image32::Set(pred, NIL);
      for(p = 0; p < n; p++)
	cost[p] = FLT_MAX;
      
      for(i=1; i<=S[0]; i++){
	cost[S[i]] = 0.0;
	Heap::Insert_MinPolicy(Q, S[i]);
      }
	
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq[i];
	      
	      if(tmp < cost[q]){
		Heap::Update_MinPolicy(Q, q, tmp);
		pred->data[q] = p;
	      }
	    }
	  }
	}
      }
      free(Dpq);
      gft::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }


    
    //------------------------------------

    //Outer Cut:
    void OIFT(sImage32 *W,
	      sAdjRel *A,
	      sImage32 *img,
	      float per,
	      int *S,
	      sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w,Wmax;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      float per_pq;
      
      value = Image32::Create(W->ncols,
			      W->nrows);
      n = label->n;
      Wmax = gft::Image32::GetMaxVal(W)*2;
      Wmax *= (1.0 + fabsf(per)/100.0);     
      Q = PQueue32::Create(Wmax+2,n,value->data);

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){

	      w = W->data[p] + W->data[q];
	      if(label->data[p] > 0)
		per_pq = per;
	      else
		per_pq = -per;

	      if(img->data[p] > img->data[q])
		w *= (1.0 + per_pq/100.0);
	      else if(img->data[p] < img->data[q])
		w *= (1.0 - per_pq/100.0);
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }



    //---------------------------------------
    
    void OIFT_in(sImageGraph *sg,
		 int *S,
		 sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }

      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){	  
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p] != 0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }
    
    
    void OIFT(sImageGraph *sg,
	      int *S,
	      sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }


    void OIFT_MaxMin(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,it = 0;
      int w;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      Image32::Set(pred, NIL);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+2+1,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = sg->Wmax*2+2; //sg->Wmax+1; //INT_MAX;
	else                    value->data[p] = 0; 
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  PQueue32::FastInsertElem(Q, S[i]);
	  //label->data[S[i]] += 2;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p]!=NIL){
	    PQueue32::FastInsertElem(Q, p);
	    //label->data[p] += 2;
	  }
	}
      }

      while(!PQueue32::IsEmpty(Q)) {
	if(it == niter)
	  break;

	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	//label->data[p] -= 2;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
		w = MAX(value->data[p], w*2);
	      }
	      else{
		w = (sg->n_link[p])[i];
		w = MAX(value->data[p], w*2+1);
	      }
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
		//label->data[q] += 2;
	      }
	    }
	  }
	}

	it++;
      }
      /*
      for(u_y = 0; u_y < label->nrows; u_y++){
	for(u_x = 0; u_x < label->ncols; u_x++)
	  printf("%d ", label->array[u_y][u_x]);
	printf("\n");
      }
      */      
      PQueue32::Destroy(&Q);
      free(i_inv);
    }

    

    void OIFT_MinMax(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *value,
		     int niter){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,it = 0;
      int w;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      Image32::Set(pred, NIL);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = 0;
	else                    value->data[p] = sg->Wmax+1; //INT_MAX;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  PQueue32::FastInsertElem(Q, S[i]);
	  label->data[S[i]] += 2;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p]!=NIL){
	    PQueue32::FastInsertElem(Q, p);
	    label->data[p] += 2;
	  }
	}
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMaxFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	label->data[p] -= 2;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w > value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
		label->data[q] += 2;
	      }
	    }
	  }
	}

	it++;
	if(it == niter)
	  break;
      }
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_TZ2Bkg(sImageGraph *sg,
		     int *S,
		     sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+3,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j] * 2;
	      }
	      else
		w = (sg->n_link[p])[i] * 2 + 1;

	      w = MAX(w, value->data[p]);
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_TZ2Obj(sImageGraph *sg,
		     int *S,
		     sImage32 *label){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax*2+3,n,value->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j] * 2 + 1;
	      }
	      else
		w = (sg->n_link[p])[i] * 2;

	      w = MAX(w, value->data[p]);
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
      free(i_inv);
    }



    

    void OIFT_TZ(sImageGraph *sg,
		 int *S,
		 sImage32 *label){
      sImage32 *lb_tzbkg, *lb_tzobj;
      int p;
      lb_tzbkg = gft::Image32::Clone(label);
      OIFT_TZ2Bkg(sg, S, lb_tzbkg);
      gft::ImageGraph::Transpose(sg);
      lb_tzobj = gft::Image32::Create(label);
      for(p = 0; p < label->n; p++){
	if(label->data[p] == NIL)
	  lb_tzobj->data[p] = NIL;
	else if(label->data[p] == 0)
	  lb_tzobj->data[p] = 1;
	else
	  lb_tzobj->data[p] = 0;
      }
      OIFT_TZ2Bkg(sg, S, lb_tzobj);
      for(p = 0; p < label->n; p++){
	if(lb_tzbkg->data[p] == 1)
	  label->data[p] = 1;
	else if(lb_tzobj->data[p] == 1)
	  label->data[p] = 0;
	else
	  label->data[p] = 2; /*Tie-zone*/
      }
      gft::ImageGraph::Transpose(sg);
      gft::Image32::Destroy(&lb_tzbkg);
      gft::Image32::Destroy(&lb_tzobj);
    }



    bool isOIFT(sImageGraph *sg,
		int *S,
		sImage32 *Slabel,
		sImage32 *label,
		sImage32 *pred,
		sImage32 *ord,
		bool complete_check){
      bool seg, forest;
      seg = isOIFT_Segmentation(sg,
				S,
				Slabel,
				label);
      if(complete_check)
	forest = isOIFT_Forest(sg,
			       S,
			       Slabel,
			       pred,
			       ord);
      else
	forest = isForest(sg,
			  pred);
      return (seg && forest);
    }

    
    
    bool isOIFT_Segmentation(sImageGraph *sg,
			     int *S,
			     sImage32 *Slabel,
			     sImage32 *label){
      //*************
      //static int i = 0;
      char filename[512];
      //*************
      gft::sImage32 *tz;
      bool flag = true;
      bool energy_test;
      int p;
      tz = gft::Image32::Clone(Slabel);
      OIFT_TZ(sg, S, tz);
      //*************
      //sprintf(filename, "tiezone.pgm");
      //Image32::Write(tz, filename);
      //*************
      for(p = 0; p < label->n; p++){
	if(label->data[p] != tz->data[p] &&
	   tz->data[p] != 2){
	  flag = false;
	  break;
	}
      }
      energy_test = (GetEnergy_Min(sg, tz, 1) == GetEnergy_Min(sg, label, 1));
      printf("within tie-zone: %d, energy test: %d\n", flag, energy_test);
      gft::Image32::Destroy(&tz);
      return (flag && energy_test);
    }


    

    bool isOIFT_Forest_tmp(sImageGraph *sg,
			   sImage32  *label,
			   sImage32  *Tpred,
			   sImage32  *value,
			   sPQueue32 *Q,
			   int *i_inv,
			   sImage32  *pred,
			   sImage32 *ord){
      struct node_oift_info { int label; int pred; int value; int color; };
      struct node_oift_info *backup;
      int u_x,u_y,v_x,v_y;
      int i,j,k,p,q,n,np,x;
      int w,bucket;
      bool flag = true;
      sAdjRel *A;
      int *F = NULL;
      A = sg->A;
      n = label->n;
      if(PQueue32::IsEmpty(Q))
	return true;
      
      bucket = Q->C.minvalue;
      while(Q->C.first[bucket] == NIL)
	bucket++;
      
      np = 0;
      p = Q->C.first[bucket];
      do{
	if(Tpred->data[p] == pred->data[p])
	  np++;
	p = Q->L.elem[p].next;
      }while(p != NIL);

      //printf("np: %d\n", np);
      
      if(np == 0)
	return false;
      
      F = (int *)calloc(np, sizeof(int));	
      np = 0;
      p = Q->C.first[bucket];
      w = value->data[p];
      do{
	if(Tpred->data[p] == pred->data[p]){
	  F[np] = p;
	  np++;
	}
	p = Q->L.elem[p].next;
      }while(p != NIL);	

      //ordenacao por insercao:
      for(i = 0; i < np-1; i++){
	// Insere F[i+1] em F[0],...,F[i].
	x = F[i+1];
	j = i;
	while( j >= 0 && ord->data[F[j]] > ord->data[x]){
	  F[j+1] = F[j];
	  j -= 1;
	}
	F[j+1] = x;
      }
      
      backup = (struct node_oift_info *)calloc(A->n,
					       sizeof(struct node_oift_info));
      
      for(k = 0; k < np; k++){
	p = F[k];
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);

	//Backup:
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    backup[i].value = value->data[q];
	    backup[i].label = label->data[q];
	    backup[i].pred  = Tpred->data[q];
	    backup[i].color = Q->L.elem[q].color;
	  }
	}

	PQueue32::FastRemoveElem(Q, p);
	flag = true;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = w;
		label->data[q] = label->data[p];
		Tpred->data[q] = p;
		PQueue32::FastInsertElem(Q, q);
	      }
	      
	      if(Tpred->data[q] != p && p == pred->data[q]){
		flag = false;
		break;
	      }
	      
	    }
	  }
	}

	if(flag)
	  flag = isOIFT_Forest_tmp(sg, label, Tpred, value,
				   Q, i_inv, pred, ord);
	
	if(flag)
	  break;

	//Restore:
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    label->data[q] = backup[i].label;
	    Tpred->data[q] = backup[i].pred;
	    if(backup[i].color == WHITE){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value->data[q] = backup[i].value;
	    }
	    else if(backup[i].color == GRAY){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value->data[q] = backup[i].value;
	      PQueue32::FastInsertElem(Q, q);
	    }
	    Q->L.elem[q].color = backup[i].color;
	  }
	}
	PQueue32::FastInsertElem(Q, p);
      }
      free(F);
      free(backup);
      
      return flag;
    }

    

    
    bool isOIFT_Forest(sImageGraph *sg,
		       int *S,
		       sImage32 *Slabel,
		       sImage32 *pred,
		       sImage32 *ord){
      sPQueue32 *Q=NULL;
      int i,p,n;
      sImage32 *value, *label, *Tpred;
      int *i_inv;
      bool flag = true;

      if(!isForest(sg, pred))
	return false;
      
      label = Image32::Clone(Slabel);
      Tpred = Image32::Create(sg->ncols,
			      sg->nrows);
      Image32::Set(Tpred, NIL);
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);

      i_inv = gft::AdjRel::InverseIndexes(sg->A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      flag = isOIFT_Forest_tmp(sg,
			       label,
			       Tpred,
			       value,
			       Q,
			       i_inv,
			       pred,
			       ord);

      Image32::Destroy(&value);
      Image32::Destroy(&label);
      Image32::Destroy(&Tpred);
      PQueue32::Destroy(&Q);
      free(i_inv);      
      return flag;
    }



    bool isForest(sImageGraph *sg,
		  sImage32 *pred){
      gft::sBMap *path,*acyclic;
      bool forest = true;
      int p,q;
      acyclic = gft::BMap::Create(pred->n);
      gft::BMap::Fill(acyclic, 0);
      path = gft::BMap::Create(pred->n);
      gft::BMap::Fill(path, 0);      
      for(p = 0; p < pred->n; p++){
	q = p;
	do{
	  if(gft::BMap::Get(path, q) == 1){
	    forest = false;
	    break;
	  }
	  else if(gft::BMap::Get(acyclic, q) == 1)
	    break;

	  gft::BMap::Set1(path, q);
	  
	  q = pred->data[q];
	}while(q != NIL);

	if(!forest) break;

	q = p;
	do{
	  if(gft::BMap::Get(acyclic, q) == 1)
	    break;
	  gft::BMap::Set0(path, q);
	  gft::BMap::Set1(acyclic, q);	  
	  q = pred->data[q];
	}while(q != NIL);
      }      
      gft::BMap::Destroy(&acyclic);
      gft::BMap::Destroy(&path);

      if(!forest)
	printf("Cycle detected.\n");
      
      return forest;
    }




    void OIFT_guided(sImageGraph *sg,
		     int *S,
		     sImage32 *label,
		     sImage32 *pred,
		     sImage32 *ord){
      sHeap32fi_lex *Q=NULL;
      int i,j,p,q,n;
      int w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      value = (float *)calloc(label->n, sizeof(float));
      n = label->ncols*label->nrows;
      Q = gft::Heap32fi_lex::Create(label->n, value, ord->data);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);

      Image32::Set(pred, NIL);
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value[p] = FLT_MAX;
	else                    value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  gft::Heap32fi_lex::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    gft::Heap32fi_lex::Insert_MinPolicy(Q, p);
      }

      while(!Heap32fi_lex::IsEmpty(Q)) {
	Heap32fi_lex::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      if(label->data[p]==0){
		j = i_inv[i];
		w = (sg->n_link[q])[j];
	      }
	      else
		w = (sg->n_link[p])[i];
	      
	      if(w < value[q]){
		//value[q] = w;
		label->data[q] = label->data[p];
		pred->data[q] = p;
		gft::Heap32fi_lex::Update_MinPolicy(Q, q, w, ord->data[q]);
	      }
	    }
	  }
	}
      }
      free(value);
      Heap32fi_lex::Destroy(&Q);
      free(i_inv);
    }


    
    
    
    void OIFT(sGraph *graph,
	      sGraph *transpose,
	      int *S,
	      int *label){
      sPQueue32 *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int w;
      int *value;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));
      Q = PQueue32::Create(Wmax+2, n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = INT_MAX;
	else              value[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);

	if(label[p]==0) g = transpose;
	else   	        g = graph;

	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->L.elem[q].color != BLACK){

	    /*
	    if(label[p]==0)
	      w = Graph::GetArcWeight(graph, q, p);
	    else
            */
	    w = g->nodes[p].Warcs[i];
	    
	    if(w < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      value[q] = w;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Q);
    }
    



    void OIFT_Heap(sImageGraph *sg,
		   int *S,
		   sImage32 *label){
      sHeap *Q=NULL;
      int i,j,p,q,n;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;

      n = label->ncols*label->nrows;
      value = gft::AllocFloatArray(n);
      Q = Heap::Create(n, value);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value[p] = FLT_MAX;
	else                    value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    Heap::Insert_MinPolicy(Q, p);
      }

      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v_x,v_y)){
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Q->color[q] != BLACK){
      
	      if(label->data[p]==0){
		j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		wi = (sg->n_link[q])[j];
	      }
	      else
		wi = (sg->n_link[p])[i];

	      if(wi == INT_MAX)
		w = FLT_MAX;
	      else
		w = (float)wi;
	      
	      if(w < value[q]){
		label->data[q] = label->data[p];
		Heap::Update_MinPolicy(Q, q, w);
	      }
	    }
	  }
	}
      }
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Q);
      free(i_inv);
    }



    void OIFT_Heap(sGraph *graph,
		   sGraph *transpose,
		   int *S,
		   int *label){
      sHeap *Q=NULL;
      sGraph *g;
      int i,j,p,q,n;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;

      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      Q = Heap::Create(n, value);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = FLT_MAX;
	else              value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    Heap::Insert_MinPolicy(Q, p);
      }

      while(!Heap::IsEmpty(Q)) {
	Heap::Remove_MinPolicy(Q, &p);

	if(label[p]==0) g = transpose;
	else   	        g = graph;
	
	for(i = 0; i < g->nodes[p].outdegree; i++){
	  q = g->nodes[p].adjList[i];

	  if(Q->color[q] != BLACK){

	    /*
	    if(label[p]==0)
	      wi = Graph::GetArcWeight(graph, q, p);
	    else
	    */
	    wi = g->nodes[p].Warcs[i];
	    
	    if(wi == INT_MAX)
	      w = FLT_MAX;
	    else
	      w = (float)wi;
	    
	    if(w < value[q]){
	      label[q] = label[p];
	      Heap::Update_MinPolicy(Q, q, w);
	    }
	  }
	}
      }
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Q);
    }
    
    

    /*
    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant,e_obj,e_bkg,e_max;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Qobj = PQueue32::Create(sg->Wmax+2,n,value->data);
      Qbkg = PQueue32::Create(sg->Wmax+2,n,value->data);
      Q	= gft::AllocIntArray(n);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value->data[S[i]] = 0;
	  if(label->data[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value->data[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)) {
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value->data[p_obj];
	e_bkg = value->data[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		
		if(label->data[p]==0){
		  j = i_inv[i];
		  w = (sg->n_link[q])[j];
		}
		else
		  w = (sg->n_link[p])[i];

                if(w < e_max){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  else if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = BLACK;
		  Qbkg->L.elem[q].color = BLACK;
		  label->data[q] = label->data[p];
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value->data[q]){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = WHITE;
		  Qbkg->L.elem[q].color = WHITE;
		  value->data[q] = w;
		  label->data[q] = label->data[p];
		  if(label->data[q] > 0)
		    PQueue32::FastInsertElem(Qobj, q);
		  else
		    PQueue32::FastInsertElem(Qbkg, q);
		}
	      }
	    }
	  }
	}
      }

      while(!PQueue32::IsEmpty(Qobj)){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){

		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		Qobj->L.elem[q].color = BLACK;
		label->data[q] = label->data[p];
		Qtop++;
		Q[Qtop] = q;
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    */

    void EOIFT(sImageGraph *sg,
	       int *S,
	       sImage32 *label,
	       int e_max){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant,e_obj,e_bkg;
      int w;
      sImage32 *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Qobj = PQueue32::Create(sg->Wmax+2,n,value->data);
      Qbkg = PQueue32::Create(sg->Wmax+2,n,value->data);
      Q	= gft::AllocIntArray(n);
      A = sg->A;

      i_inv = gft::AdjRel::InverseIndexes(A);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value->data[S[i]] = 0;
	  if(label->data[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value->data[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value->data[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value->data[p] = INT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      //while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)){
      while(Qobj->nadded != 0 && Qbkg->nadded != 0){
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value->data[p_obj];
	e_bkg = value->data[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	lp = label->data[p];

	goto label04;
	
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;

	label04:

	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  w = (sg->n_link[q])[j];
		}
		else
		  w = (sg->n_link[p])[i];

                if(w < e_max){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  else if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = BLACK;
		  Qbkg->L.elem[q].color = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value->data[q]){
		  if(Qobj->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qobj, q);
		  if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		  Qobj->L.elem[q].color = WHITE;
		  Qbkg->L.elem[q].color = WHITE;
		  value->data[q] = w;
		  label->data[q] = lp;
		  if(lp > 0)
		    PQueue32::FastInsertElem(Qobj, q);
		  else
		    PQueue32::FastInsertElem(Qbkg, q);
		}
	      }
	    }
	  }
	}
      }

      //while(!PQueue32::IsEmpty(Qobj)){
      while(Qobj->nadded != 0){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->L.elem[q].color != BLACK){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		Qobj->L.elem[q].color = BLACK;
		label->data[q] = label->data[p];
		Qtop++;
		Q[Qtop] = q;
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    


    //------------

    void EOIFT(sGraph *graph,
	       sGraph *transpose,
	       int *S,
	       int *label,
	       int e_max){
      sPQueue32 *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant,e_obj,e_bkg;
      int w;
      int *value;
      int *Q=NULL;
      int Qtop = -1;
      int Wmax;
      Wmax = Graph::GetMaximumArc(graph);
      n = graph->nnodes;
      value = (int *)malloc(n*sizeof(int));

      Qobj = PQueue32::Create(Wmax+2, n, value);
      Qbkg = PQueue32::Create(Wmax+2, n, value);
      Q	= gft::AllocIntArray(n);

      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0;
	  if(label[S[i]] == 0)
	    PQueue32::FastInsertElem(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    PQueue32::FastInsertElem(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = INT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0;
	    PQueue32::FastInsertElem(Qbkg, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0;
	    PQueue32::FastInsertElem(Qobj, p);
	  }
	  else{
	    value[p] = INT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      //while(!PQueue32::IsEmpty(Qobj) && !PQueue32::IsEmpty(Qbkg)){
      while(Qobj->nadded != 0 && Qbkg->nadded != 0){
	p_obj = PQueue32::FastGetMinFIFO(Qobj);
	p_bkg = PQueue32::FastGetMinFIFO(Qbkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  p = p_obj;
	  PQueue32::FastRemoveElem(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  p = p_bkg;
	  PQueue32::FastRemoveElem(Qbkg, p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    p = p_obj;
	    PQueue32::FastRemoveElem(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    PQueue32::FastRemoveElem(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}
	
	Qobj->L.elem[p].color = BLACK;
	Qbkg->L.elem[p].color = BLACK;

	lp = label[p];

	if(lp==0) g = transpose;
	else      g = graph;
	
	goto label05;

	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;

	label05:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(Qobj->L.elem[q].color != BLACK){
	      w = g->nodes[p].Warcs[i];
	      
	      if(w < e_max){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		if(Qbkg->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Qbkg, q);
		Qobj->L.elem[q].color = BLACK;
		Qbkg->L.elem[q].color = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		if(Qobj->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qobj, q);
		if(Qbkg->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Qbkg, q);
		Qobj->L.elem[q].color = WHITE;
		Qbkg->L.elem[q].color = WHITE;
		value[q] = w;
		label[q] = lp;
		if(lp > 0)
		  PQueue32::FastInsertElem(Qobj, q);
		else
		  PQueue32::FastInsertElem(Qbkg, q);		  
	      }
	      
	    }
	  }
	}
      }
      
      //while(!PQueue32::IsEmpty(Qobj)){
      while(Qobj->nadded != 0){
	p = PQueue32::FastRemoveMinFIFO(Qobj);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    
	    if(Qobj->L.elem[q].color != BLACK){
	      
	      if(Qobj->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Qobj, q);
	      Qobj->L.elem[q].color = BLACK;
	      label[q] = label[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      free(value);
      PQueue32::Destroy(&Qobj);
      PQueue32::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
    }
    

    //----------------------------------------
    /*
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant;
      float e_obj,e_bkg,e_max;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      //----------
      n = label->n;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= gft::AllocIntArray(n);
      A = sg->A;
      
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(!Heap::IsEmpty(Qobj) && !Heap::IsEmpty(Qbkg)) {
	Heap::Get_MinPolicy(Qobj, &p_obj);
	Heap::Get_MinPolicy(Qbkg, &p_bkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	  Heap::Delete_MinPolicy(Qobj, p);
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	  Heap::Delete_MinPolicy(Qbkg, p);
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0){
	    p = p_obj;
	    Heap::Delete_MinPolicy(Qobj, p);
	  }
          else{
	    p = p_bkg;
	    Heap::Delete_MinPolicy(Qbkg, p);
	  }
          l_ant = 1 - l_ant;
	}

	//-----------------
	//printf("e_obj: %5d, e_bkg: %5d, ", ROUND(e_obj), ROUND(e_bkg));
	//printf("x: %4d, y: %4d\n", p%label->ncols, p/label->ncols);
	//-----------------
	
	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	//Stack::Push(Q, p);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  //p = Stack::Pop(Q);
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols; //PixelX(label, p);
	  u_y = p/label->ncols; //PixelY(label, p);
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){
		
		if(label->data[p]==0){
		  j = i_inv[i]; //j = ImageGraph::get_edge_index(q, p, sg);
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX)
		  w = FLT_MAX;
		else
		  w = (float)wi;
		
                if(w < e_max){
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  else if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Qobj->color[q] = BLACK;
		  Qbkg->color[q] = BLACK;
		  label->data[q] = label->data[p];
		  //Stack::Push(Q, q);
		  Qtop++;
		  Q[Qtop] = q;
		  //---------
		  //NbyQ++;
		}
		else if(w < value[q]){
		  label->data[q] = label->data[p];
		  if(label->data[q] > 0){
		    if(Qbkg->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qbkg, q);
		    Heap::Update_MinPolicy(Qobj, q, w);
		  }
		  else{
		    if(Qobj->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qobj, q);
		    Heap::Update_MinPolicy(Qbkg, q, w);
		  }
		}
	      }
	    }
	  }
	}
      }

      while(!Heap::IsEmpty(Qobj)){
	Heap::Remove_MinPolicy(Qobj, &p);
	//Stack::Push(Q, p);
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  //p = Stack::Pop(Q);
	  p = Q[Qtop];
	  Qtop--;
	  u_x = p%label->ncols; //PixelX(label, p);
	  u_y = p/label->ncols; //PixelY(label, p);
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(Image32::IsValidPixel(label,v_x,v_y)){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){

		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		Qobj->color[q] = BLACK;
		label->data[q] = label->data[p];
		//Stack::Push(Q, q);
		Qtop++;
		Q[Qtop] = q;
		//--------
		//NbyQ++;
	      }
	    }
	  }
	}
      }

      //printf("NbyQ: %d -> %f\n",NbyQ, (float)NbyQ/(float)n);
      
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    */
    //------------------------------------

    /*
    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      //----------
      n = label->n;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= gft::AllocIntArray(n);

      A = sg->A;
      
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label->data[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(Qobj->last > 0 && Qbkg->last > 0){
      //---------------------
	p_obj = Qobj->pixel[1];
	p_bkg = Qbkg->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  Heap::Remove_MinPolicy(Qobj, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  Heap::Remove_MinPolicy(Qbkg, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0){
	    Heap::Remove_MinPolicy(Qobj, &p);
	  }
          else{
	    Heap::Remove_MinPolicy(Qbkg, &p);
	  }
          l_ant = 1 - l_ant;
	}

	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	lp = label->data[p];
	
	goto label01;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label01:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(Qobj->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  else if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Qobj->color[q] = BLACK;
		  Qbkg->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value[q]){
		  label->data[q] = lp;
		  if(lp > 0){
		    if(Qbkg->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qbkg, q);
		    Heap::Update_MinPolicy(Qobj, q, w);
		  }
		  else{
		    if(Qobj->color[q] == GRAY)
		      Heap::Delete_MinPolicy(Qobj, q);
		    Heap::Update_MinPolicy(Qbkg, q, w);
		  }
		}
	      }
	    }
	  }
	}
      }

      for(i = 1; i <= Qobj->last; i++){
	p = Qobj->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	Qobj->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(Qobj->color[q] != BLACK){
	      Qobj->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    */
    //-----------------------------------------------------


    void EOIFT_Heap(sImageGraph *sg,
		    int *S,
		    sImage32 *label,
		    float e_max){
      sHeapPair *QP=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      n = label->n;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      A = sg->A;
      
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label->data[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label->data[p] != NIL){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0)
	    HeapPair::Remove_MinPolicy_1(QP, &p);
          else
	    HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	QP->color[p] = BLACK;

	lp = label->data[p];
	
	goto label03;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label03:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(QP->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(QP->color[q] == GRAY)
		    HeapPair::Delete_MinPolicy(QP, q);
		  QP->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < value[q]){
		  label->data[q] = lp;
		  if(lp > 0)
		    HeapPair::Update_MinPolicy_1(QP, q, w);
		  else
		    HeapPair::Update_MinPolicy_0(QP, q, w);
		}
	      }
	    }
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(QP->color[q] != BLACK){
	      QP->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
      free(i_inv);
    }
    

    //-----------------------------------------------------
    
    void EOIFT_Heap_2(sImageGraph *sg,
		      int *S,
		      sImage32 *label,
		      float e_max){
      sHeapPair *QP=NULL;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sAdjRel *A;
      int *i_inv;
      int *Q=NULL;
      int Qtop = -1;
      int *T=NULL;
      float *Tv = NULL;
      int Ttop = -1;
      n = label->n;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      T	= gft::AllocIntArray(n);
      Tv = gft::AllocFloatArray(n);

      A = sg->A;
      i_inv = gft::AdjRel::InverseIndexes(A);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  Tv[S[i]] = value[S[i]] = 0.0;
	  if(label->data[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label->data[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label->data[p]==NIL){
	    Tv[p] = value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	  else
	    Tv[p] = value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label->data[p] == 0){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label->data[p] != NIL){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    Tv[p] = value[p] = FLT_MAX;
	    label->data[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	  //------------------
	  //p_obj = QP->pixel[QP->n];
	  //while(QP->last_1 <= QP->n && value[p_obj] < e_max){
	  //  HeapPair::Remove_MinPolicy_1(QP, &p_obj);
	  //  Qtop++;
	  //  Q[Qtop] = p_obj;
	  //  p_obj = QP->pixel[QP->n];
	  //}
	  //------------------
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	  //------------------
	  //p_bkg = QP->pixel[1];
	  //while(QP->last_0 > 0 && value[p_bkg] < e_max){
	  //  HeapPair::Remove_MinPolicy_0(QP, &p_bkg);
	  //  Qtop++;
	  //  Q[Qtop] = p_bkg;
	  //  p_bkg = QP->pixel[1];
	  //}
	  //------------------
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else           HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	  //------------------
	  //if(e_obj < e_max){
	  //  p_obj = QP->pixel[QP->n];
	  //  while(QP->last_1 <= QP->n && value[p_obj] < e_max){
	  //    HeapPair::Remove_MinPolicy_1(QP, &p_obj);
	  //    Qtop++;
	  //    Q[Qtop] = p_obj;
	  //    p_obj = QP->pixel[QP->n];
	  //  }
	  //  p_bkg = QP->pixel[1];
	  //  while(QP->last_0 > 0 && value[p_bkg] < e_max){
	  //    HeapPair::Remove_MinPolicy_0(QP, &p_bkg);
	  //    Qtop++;
	  //    Q[Qtop] = p_bkg;
	  //    p_bkg = QP->pixel[1];
	  //  }
	  //}
	  //else e_max = e_obj;
	  //------------------
	}

	//QP->color[p] = BLACK;

	lp = label->data[p];
	
	goto label01;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label01:
	  
	  u_x = p%label->ncols;
	  u_y = p/label->ncols;
	  for(i=1; i<A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if(v_x >= 0 && v_x < label->ncols &&
	       v_y >= 0 && v_y < label->nrows){
	      q = v_x + label->ncols*v_y;
	      if(QP->color[q] != BLACK){
		
		if(lp == 0){
		  j = i_inv[i];
		  wi = (sg->n_link[q])[j];
		}
		else
		  wi = (sg->n_link[p])[i];

		if(wi == INT_MAX) continue;

		w = (float)wi;
		
                if(w < e_max){
		  if(QP->color[q] == GRAY)
		    HeapPair::Delete_MinPolicy(QP, q);
		  QP->color[q] = BLACK;
		  label->data[q] = lp;
		  Qtop++;
		  Q[Qtop] = q;
		}
		else if(w < Tv[q]){ //if(w < value[q]){
		  if(Tv[q] >= value[q]){
		    label->data[q] = lp;
		    Ttop++;
		    T[Ttop] = q;
		  }
		  Tv[q] = w;
		}
	      }
	    }
	  }
	}
	while(Ttop > -1){
	  p = T[Ttop];
	  w = Tv[p];
	  Ttop--;
	  if(QP->color[p] != BLACK){
	    if(label->data[p] > 0)
	      HeapPair::Update_MinPolicy_1(QP, p, w);
	    else
	      HeapPair::Update_MinPolicy_0(QP, p, w);
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }
      
      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;
	u_x = p%label->ncols;
	u_y = p/label->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < label->ncols &&
	     v_y >= 0 && v_y < label->nrows){
	    q = v_x + label->ncols*v_y;
	    if(QP->color[q] != BLACK){
	      QP->color[q] = BLACK;
	      label->data[q] = label->data[p];
	      Qtop++;
	      Q[Qtop] = q;
	    }
	  }
	}
      }
      
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
      gft::FreeIntArray(&T);
      gft::FreeFloatArray(&Tv);
      free(i_inv);
    }


    //------------------------------------
    /*
    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n;
      int l_ant;
      float e_obj,e_bkg,e_max;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      sStack *Q=NULL;
      
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= Stack::Create(n);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) value[p] = FLT_MAX;
	else              value[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  if(label[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p] == 0)
	    Heap::Insert_MinPolicy(Qbkg, p);
	  else if(label[p] != NIL)
	    Heap::Insert_MinPolicy(Qobj, p);
      }

      l_ant = 0;
      while(!Heap::IsEmpty(Qobj) && !Heap::IsEmpty(Qbkg)) {
	Heap::Get_MinPolicy(Qobj, &p_obj);
	Heap::Get_MinPolicy(Qbkg, &p_bkg);

	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = e_bkg;
	  p = p_obj;
	}
	else if(e_obj > e_bkg){
	  e_max = e_obj;
	  p = p_bkg;
	}
	else{
	  e_max = e_obj;
          if(l_ant == 0)
	    p = p_obj;
          else
	    p = p_bkg;
          l_ant = 1 - l_ant;
	}

	if(Qobj->color[p] == GRAY)
	  Heap::Delete_MinPolicy(Qobj, p);
	if(Qbkg->color[p] == GRAY)
	  Heap::Delete_MinPolicy(Qbkg, p);
	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	Stack::Push(Q, p);
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);
	  
	  if(label[p]==0) g = transpose;
	  else   	  g = graph;
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];

	    if(Qobj->color[q] != BLACK){

	      
	      //if(label[p]==0)
		//wi = Graph::GetArcWeight(graph, q, p);
	      //else
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX)
		w = FLT_MAX;
	      else
		w = (float)wi;

	      if(w < e_max){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		Qobj->color[q] = BLACK;
		Qbkg->color[q] = BLACK;
		
		label[q] = label[p];
		Stack::Push(Q, q);
	      }
	      else if(w < value[q]){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		value[q] = w;
		label[q] = label[p];
		if(label[q] > 0)
		  Heap::Insert_MinPolicy(Qobj, q);
		else
		  Heap::Insert_MinPolicy(Qbkg, q);		  
	      }
	      
	    }
	  }

	}
	
      }

      while(!Heap::IsEmpty(Qobj)){
	Heap::Remove_MinPolicy(Qobj, &p);
	Stack::Push(Q, p);
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);

	  for(i = 0; i < graph->nodes[p].outdegree; i++){
	    q = graph->nodes[p].adjList[i];
	    
	    if(Qobj->color[q] != BLACK){

	      if(Qobj->color[q] == GRAY)
		Heap::Delete_MinPolicy(Qobj, q);
	      Qobj->color[q] = BLACK;
	      label[q] = label[p];
	      Stack::Push(Q, q);
	    }
	  }
	}
      }

      while(!Heap::IsEmpty(Qbkg)){
	Heap::Remove_MinPolicy(Qbkg, &p);
	Stack::Push(Q, p);
	Qobj->color[p] = BLACK;
	while(!Stack::IsEmpty(Q)){
	  p = Stack::Pop(Q);

	  for(i = 0; i < transpose->nodes[p].outdegree; i++){
	    q = transpose->nodes[p].adjList[i];
	    
	    if(Qobj->color[q] != BLACK){
	      
	      if(Qbkg->color[q] == GRAY)
		Heap::Delete_MinPolicy(Qbkg, q);
	      Qobj->color[q] = BLACK;
	      label[q] = label[p];
	      Stack::Push(Q, q);
	    }
	  }
	}
      }

      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      Stack::Destroy(&Q);
    }
    */



    /*    
    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max){
      sHeap *Qobj=NULL, *Qbkg=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int u_x,u_y,v_x,v_y;
      int *Q=NULL;
      int Qtop = -1;
      
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      Qobj = Heap::Create(n, value);
      Qbkg = Heap::Create(n, value);
      Q	= gft::AllocIntArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    Heap::Insert_MinPolicy(Qbkg, S[i]);
	  else if(label[S[i]] != NIL)
	    Heap::Insert_MinPolicy(Qobj, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qbkg, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0.0;
	    Heap::Insert_MinPolicy(Qobj, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(Qobj->last > 0 && Qbkg->last > 0){
	p_obj = Qobj->pixel[1];
	p_bkg = Qbkg->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  Heap::Remove_MinPolicy(Qobj, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  Heap::Remove_MinPolicy(Qbkg, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0)
	    Heap::Remove_MinPolicy(Qobj, &p);
          else
	    Heap::Remove_MinPolicy(Qbkg, &p);
          l_ant = 1 - l_ant;
	}

	Qobj->color[p] = BLACK;
	Qbkg->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label02;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];

	    if(Qobj->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(Qobj->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qobj, q);
		else if(Qbkg->color[q] == GRAY)
		  Heap::Delete_MinPolicy(Qbkg, q);
		Qobj->color[q] = BLACK;
		Qbkg->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		label[q] = lp;
		if(lp > 0){
		  if(Qbkg->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qbkg, q);
		  Heap::Update_MinPolicy(Qobj, q, w);
		}
		else{
		  if(Qobj->color[q] == GRAY)
		    Heap::Delete_MinPolicy(Qobj, q);
		  Heap::Update_MinPolicy(Qbkg, q, w);
		}
	      }
	      
	    }
	  }

	}
	
      }

      for(i = 1; i <= Qobj->last; i++){
	p = Qobj->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	Qobj->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Qobj->color[q] != BLACK){
	    Qobj->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      gft::FreeFloatArray(&value);
      Heap::Destroy(&Qobj);
      Heap::Destroy(&Qbkg);
      gft::FreeIntArray(&Q);
    }
    */



    void EOIFT_Heap(sGraph *graph,
		    sGraph *transpose,
		    int *S,
		    int *label,
		    float e_max){
      sHeapPair *QP=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int *Q=NULL;
      int Qtop = -1;
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label[p] != NIL){
	    value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else   	 HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	QP->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label06;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label06:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(QP->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(QP->color[q] == GRAY)
		  HeapPair::Delete_MinPolicy(QP, q);
		QP->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < value[q]){
		label[q] = lp;
		if(lp > 0)
		  HeapPair::Update_MinPolicy_1(QP, q, w);
		else
		  HeapPair::Update_MinPolicy_0(QP, q, w);
	      }
	    }
	  }
	}
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(QP->color[q] != BLACK){
	    QP->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
    }
    

    
    void EOIFT_Heap_2(sGraph *graph,
		      sGraph *transpose,
		      int *S,
		      int *label,
		      float e_max){
      sHeapPair *QP=NULL;
      sGraph *g;
      int i,j,p,p_obj,p_bkg,q,n,lp;
      int l_ant;
      float e_obj,e_bkg;
      int wi;
      float w;
      float *value;
      int *Q=NULL;
      int Qtop = -1;
      int *T=NULL;
      float *Tv = NULL;
      int Ttop = -1;
      n = graph->nnodes;
      value = gft::AllocFloatArray(n);
      QP = HeapPair::Create(n, value);
      Q	= gft::AllocIntArray(n);
      T	= gft::AllocIntArray(n);
      Tv = gft::AllocFloatArray(n);
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  Tv[S[i]] = value[S[i]] = 0.0;
	  if(label[S[i]] == 0)
	    HeapPair::Insert_MinPolicy_0(QP, S[i]);
	  else if(label[S[i]] != NIL)
	    HeapPair::Insert_MinPolicy_1(QP, S[i]);
	}
	for(p=0; p<n; p++){
	  if(label[p]==NIL){
	    Tv[p] = value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	  else
	    Tv[p] = value[p] = 0.0;
	}
      }
      else{
	for(p=0; p<n; p++){
	  if(label[p] == 0){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_0(QP, p);
	  }
	  else if(label[p] != NIL){
	    Tv[p] = value[p] = 0.0;
	    HeapPair::Insert_MinPolicy_1(QP, p);
	  }
	  else{
	    Tv[p] = value[p] = FLT_MAX;
	    label[p] = 0;
	  }
	}
      }

      l_ant = 0;
      while(QP->last_0 > 0 && QP->last_1 <= QP->n){
	p_obj = QP->pixel[QP->n];
	p_bkg = QP->pixel[1];
	
	e_obj = value[p_obj];
	e_bkg = value[p_bkg];
	if(e_obj < e_bkg){
	  e_max = MAX(e_max, e_bkg);
	  HeapPair::Remove_MinPolicy_1(QP, &p);
	}
	else if(e_obj > e_bkg){
	  e_max = MAX(e_max, e_obj);
	  HeapPair::Remove_MinPolicy_0(QP, &p);
	}
	else{
	  e_max = MAX(e_max, e_obj);
          if(l_ant == 0) HeapPair::Remove_MinPolicy_1(QP, &p);
          else   	 HeapPair::Remove_MinPolicy_0(QP, &p);
          l_ant = 1 - l_ant;
	}

	//QP->color[p] = BLACK;

	lp = label[p];
	
	if(lp==0) g = transpose;
	else   	  g = graph;
	
	goto label02;
	  
	Qtop++;
	Q[Qtop] = p;
	while(Qtop > -1){
	  p = Q[Qtop];
	  Qtop--;
	  
	label02:
	  
	  for(i = 0; i < g->nodes[p].outdegree; i++){
	    q = g->nodes[p].adjList[i];
	    if(QP->color[q] != BLACK){
	      
	      wi = g->nodes[p].Warcs[i];

	      if(wi == INT_MAX) continue;

	      w = (float)wi;

	      if(w < e_max){
		if(QP->color[q] == GRAY)
		  HeapPair::Delete_MinPolicy(QP, q);
		QP->color[q] = BLACK;
		label[q] = lp;
		Qtop++;
		Q[Qtop] = q;
	      }
	      else if(w < Tv[q]){
		if(Tv[q] >= value[q]){
		  label[q] = lp;
		  Ttop++;
		  T[Ttop] = q;
		}
		Tv[q] = w;
	      }
	    }
	  }
	}
	while(Ttop > -1){
	  p = T[Ttop];
	  w = Tv[p];
	  Ttop--;
	  if(QP->color[p] != BLACK){
	    if(label[p] > 0)
	      HeapPair::Update_MinPolicy_1(QP, p, w);
	    else
	      HeapPair::Update_MinPolicy_0(QP, p, w);
	  }
	}	
      }

      for(i = QP->last_1; i <= QP->n; i++){
	p = QP->pixel[i];
	Qtop++;
	Q[Qtop] = p;
	QP->color[p] = BLACK;
      }

      while(Qtop > -1){
	p = Q[Qtop];
	Qtop--;      
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(QP->color[q] != BLACK){
	    QP->color[q] = BLACK;
	    label[q] = label[p];
	    Qtop++;
	    Q[Qtop] = q;
	  }
	}
      }
	
      gft::FreeFloatArray(&value);
      HeapPair::Destroy(&QP);
      gft::FreeIntArray(&Q);
      gft::FreeIntArray(&T);
      gft::FreeFloatArray(&Tv);
    }
    
    
    
    
    void IFT_fmax_Heap(sGraph *graph,
		       int *S,
		       int *label,
		       float *cost){
      sHeap *Q;
      float tmp, w;
      int n,p,q,i;
      n = graph->nnodes;
      Q = gft::Heap::Create(n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = FLT_MAX;
	else              cost[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  gft::Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    gft::Heap::Insert_MinPolicy(Q, p);
      }
      
      while(!gft::Heap::IsEmpty(Q)){
	gft::Heap::Remove_MinPolicy(Q, &p);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = MAX(cost[p], w);

	    if(tmp < cost[q]){
	      gft::Heap::Update_MinPolicy(Q, q, tmp);
	      label[q] = label[p];
	    }
	  }
	}
      }
      gft::Heap::Destroy(&Q);
    }




    void IFT_fmax(sGraph *graph,
		  int *S,
		  int *label,
		  int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = MAX(cost[p], w);

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }



    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);

      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }
    
    

    
    void IFT_fw_Heap(sGraph *graph,
		     int *S,
		     int *label,
		     float *cost,
		     int *pred){
      sHeap *Q;
      float tmp, w;
      int n,p,q,i;
      n = graph->nnodes;
      Q = gft::Heap::Create(n, cost);

      for(p = 0; p < n; p++){
	pred[p] = NIL;
	if(label[p]==NIL) cost[p] = FLT_MAX;
	else              cost[p] = 0.0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  gft::Heap::Insert_MinPolicy(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    gft::Heap::Insert_MinPolicy(Q, p);
      }
      
      while(!gft::Heap::IsEmpty(Q)){
	gft::Heap::Remove_MinPolicy(Q, &p);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->color[q] != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      gft::Heap::Update_MinPolicy(Q, q, tmp);
	      label[q] = label[p];
	      pred[q] = p;
	    }
	  }
	}
      }
      gft::Heap::Destroy(&Q);
    }
   


    void IFT_fw(sGraph *graph,
		int *S,
		int *label,
		int *cost,
		int *pred){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      n = graph->nnodes;
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create(Wmax+2, n, cost);
      
      for(p = 0; p < n; p++){
	pred[p] = NIL;
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graph->nodes[p].Warcs[i];
	    tmp = w;

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      pred[q] = p;
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
    }
   
    
    //---------------------------------------
    // Convex IFT:

    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_conquer_path(int p,
			 sImageGraph *sg,
			 sImage32 *P_sum, 
			 sImage32 *V,
			 sPQueue32 *Q,
			 sImage32 *label){
      int i,q,edge;
      Pixel u,v;
      sAdjRel *A;
      
      A = sg->A;
      do{
	if(Q->L.elem[p].color == GRAY)
	  PQueue32::FastRemoveElem(Q, p);
	Q->L.elem[p].color = BLACK;
	
	label->data[p] = 1;
	
	u.x = p%label->ncols;
	u.y = p/label->ncols;
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      
	      if(edge < V->data[q] && q != P_sum->data[p]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = edge;
		label->data[q] = 1; //label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
	p = P_sum->data[p];
	if(p == NIL) break;
      }while(Q->L.elem[p].color != BLACK);
    }

    
    /* P_sum = Predecessor map obtained by the IFT fsum.*/
    void SC_prune_tree(int p,
		       sImageGraph *sg,
		       sImage32 *P_sum, 
		       sImage32 *V,
		       sPQueue32 *Q,
		       sQueue *Qfifo,
		       sImage32 *label){
      Pixel u,v;
      int i,q,edge;
      sAdjRel *A = sg->A;
      
      if(Q->L.elem[p].color == GRAY)
	PQueue32::FastRemoveElem(Q, p);
      Q->L.elem[p].color = BLACK;
      
      label->data[p] = 0;
      
      Queue::Push(Qfifo, p);
      
      //printf("Prune tree\n");
      
      while(!Queue::IsEmpty(Qfifo)){
	p = Queue::Pop(Qfifo);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label,v.x,v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(P_sum->data[q] == p){
		Queue::Push(Qfifo, q);
		
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		Q->L.elem[q].color = BLACK;
		
		label->data[q] = 0;
	      }
	      else{
		edge = (sg->n_link[p])[i];
		if(edge < V->data[q]){
		  if(Q->L.elem[q].color == GRAY)
		    PQueue32::FastRemoveElem(Q, q);
		  V->data[q] = edge;
		  label->data[q] = 0;
		  PQueue32::FastInsertElem(Q, q);
		}
	      }
	      
	    }
	  }
	}
	
      }
    }

    
    void SC_IFT(sImageGraph *sg,
		int *S,
		sImage32 *label,
		sImage32 *P_sum){
      sPQueue32 *Q=NULL;
      sQueue *Qfifo=NULL;
      sImage32 *V;
      int p,n,i;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      Qfifo = Queue::Create(n);
      
      for(p = 0; p < n; p++){
	if(label->data[p]==NIL) V->data[p] = INT_MAX;
	else                    V->data[p] = 0;
      }

      for(i = 1; i <= S[0]; i++)
	PQueue32::FastInsertElem(Q, S[i]);
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	if(label->data[p] > 0)
	  SC_conquer_path(p, sg, P_sum, V, Q, label);
	else if(label->data[p] == 0)
	  SC_prune_tree(p, sg, P_sum, V, Q, Qfifo, label);
      }
      
      Image32::Destroy(&V);
      Queue::Destroy(&Qfifo);
      PQueue32::Destroy(&Q);
    }


    //----------------------------------------

    sImage32 *Cost_fmin(sImageGraph *sg,
			int *S, int lb,
			sImage32 *label){
      sPQueue32 *Q = NULL; 
      sImage32 *V;
      int i,p,q,n, edge,tmp;
      Pixel u,v;
      sAdjRel *A;
      
      n = sg->ncols*sg->nrows;
      V = Image32::Create(sg->ncols, sg->nrows);
      Q = PQueue32::Create(sg->Wmax+2, n, V->data);
      A = sg->A;
      
      for(p=0; p<n; p++){
	if(label->data[p] == lb) V->data[p] = sg->Wmax+1;
	else                     V->data[p] = INT_MIN;
      }

      if(S != NULL){
	for(i = 1; i <= S[0]; i++){
	  if(label->data[S[i]] == lb)
	    PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p] == lb)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMaxFIFO(Q);
	u.x = p%label->ncols; 
	u.y = p/label->ncols; 
	
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  //if(Image32::IsValidPixel(label, v.x, v.y)){
	  if(v.x >= 0 && v.x < label->ncols &&
	     v.y >= 0 && v.y < label->nrows){
	    q = v.x + label->ncols*v.y;
	    if(Q->L.elem[q].color != BLACK){
	      edge = (sg->n_link[p])[i];          
	      tmp  = MIN(V->data[p], edge);
	      if(tmp > V->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		V->data[q] = tmp; //mapa de custos
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      return V;
    }

    //----------------------------------------
    

    int *GetSeedsByLabel(int *S,
			 sImage32 *label,
			 int lb){
      int *Slb = NULL;
      int l, s, i, j, size = 0;
      for(i = 1; i <= S[0]; i++){
	s = S[i];
	l = label->data[s];
	if(l == lb) size++;
      }
      Slb = gft::AllocIntArray(size+1);
      Slb[0] = size;
      i = 1;
      for(j = 1; j <= S[0]; j++){
	s = S[j];
	l = label->data[s];
	if(l == lb){
	  Slb[i] = s;
	  i++;
	}
      }
      return Slb;
    }
    
    
    int *GetAllInternalSeedsByLabel(int *S,
				    sImage32 *label,
				    int lb,
				    int *hr,
				    int nlayers){
      int *Sall, *Sl, *Stmp;
      int l,size,i,j;
      Sall = GetSeedsByLabel(S, label, lb);
      for(l=0; l < nlayers; l++){
	if(hr[l] == lb-1){
	  //Sl = GetSeedsByLabel(S, label, l+1);
	  Sl = GetAllInternalSeedsByLabel(S, label, l+1,
					  hr, nlayers);
	  size = Sl[0] + Sall[0];
	  Stmp = gft::AllocIntArray(size+1);
	  Stmp[0] = size;
	  j =  1;
	  for(i = 1; i <= Sall[0]; i++){
	    Stmp[j] = Sall[i];
	    j++;
	  }
	  for(i = 1; i <= Sl[0]; i++){
	    Stmp[j] = Sl[i];
	    j++;
	  }
	  gft::FreeIntArray(&Sall);
	  gft::FreeIntArray(&Sl);
	  Sall = Stmp;
	}
      }
      return Sall;
    }
    

    sImageGraph *GetPolarityGraph(sImageGraph *graph,
				  sCImage *cimg,
				  sImage32 *img,
				  int pol){
      sImageGraph *sg_oriented;
      sg_oriented = gft::ImageGraph::Clone(graph);
      if(cimg != NULL){
	gft::sImage32 *lumi;
	lumi = gft::CImage::Luminosity(cimg);
	gft::ImageGraph::Orient2Digraph(sg_oriented, lumi, pol);
	gft::Image32::Destroy(&lumi);
      }
      else
	gft::ImageGraph::Orient2Digraph(sg_oriented, img, pol);
      return sg_oriented;
    }
    

    
    void HL_OIFT(sImageGraph *graph,
		 sCImage *cimg,
		 sImage32 *img,
		 float radius,
		 char *hierarchy,
		 int *S,
		 sImage32 *label){
      int p;
      //-----------------------
      gft::sLayeredGraph *lg;
      gft::sImageGraph *sg, *sg_oriented;
      //------------------------
      FILE *fpHierarchy;
      int nlayers, typeL, nRelations, relationT, layer1, layer2;
      int i, j, n, nSeeds, k, s, Wmax, pol;
      int x, y, lb, ncols, nrows;
      int *lg_label;
      
      ncols = label->ncols;
      nrows = label->nrows;
      gft::ImageGraph::ChangeType(graph, DISSIMILARITY);
      
      fpHierarchy = fopen((char *)hierarchy, "r");
      if (!fpHierarchy){
	gft::Error((char *)MSG1,(char *)"Error opening file of Hierarchy Relations!");
      }
      
      /*------------- Set number of layers -------------*/
      (void)fscanf(fpHierarchy,"%d\n",&nlayers);
      //printf("NLayers =  %d\n",nlayers);
      
      /*-------------GIVEN BY USER: TYPE OF LAYER -------------*/
      /*typeLayer:
	0 = normal,
	1 = GSC,
      */
      int *typeLayer = gft::AllocIntArray(nlayers);
      int *polLayer = gft::AllocIntArray(nlayers);
      
      for(i = 0; i < nlayers; i++) {
	(void)fscanf(fpHierarchy,"%d %d\n",&typeL, &pol);
	typeLayer[i] = typeL;
	polLayer[i] = pol;
      }
      
      /*------------- Creating HIERARCHY -------------*/
      /*Relations is given by inclusion(=1) or exclusion(=2) binary relation between: (son,father) and (brothers) that defines a hierarchy, in a .txt*/
      /*Example: 1|0|1 = inclusion of layer 0 in layer 1 
	or 2|0|1 = exclusion between layer 0 and layer 1*/
      
      /*Initialization of Hierarchy*/
      int *hr = gft::AllocIntArray(nlayers);
      for(i = 0; i < nlayers; i++){
	hr[i]= -1; 
      }
      
      /*Reading a .txt fill, and set in "hr" with only inclusion cases*/
      // hr[0] = 1;  // i(0,1) = 1 0 1  // means: "1 is father of 0" or "0 is included in 1"
      (void)fscanf(fpHierarchy,"%d\n",&nRelations);
      for(i = 0; i < nRelations; i++) {
	(void)fscanf(fpHierarchy,"%d %d %d\n",&relationT,&layer1,&layer2);
	if(relationT == 1){ // it's inclusion relation
	  hr[layer1] = layer2;
	}
      }
      
      fclose(fpHierarchy);
      
      /* ------------- Create label image -------------*/
      n = ncols*nrows*nlayers;
      lg_label = gft::AllocIntArray(n);
      /*Set NIL (-1) in all matrix positions*/
      for(i = 0; i < n; i++){
	lg_label[i] = NIL;
      }
      
      nSeeds = S[0];
      int* Seeds = gft::AllocIntArray(nSeeds*nlayers + 1);
      
      k = 1;
      for(i = 1; i <= S[0]; i++){
	s = S[i];
	x = s%ncols;
	y = s/ncols;
	lb = label->data[s];
	
	if(lb != 0){
	  p = x + y*ncols + nrows*ncols*(lb-1); 
	  lg_label[p] = 1;
	  Seeds[k] = p;
	  k++;
	  /*Essa parte faltou no algoritmo do artigo:*/
	  for(j = 0; j < nlayers ; j++){
	    /*
	    if(j != lb-1 && hr[lb-1] == j){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 1;
	      Seeds[k] = p;
	      k++;
	    }
	    else if(j != lb-1){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 0;
	      Seeds[k] = p;
	      k++;
	    }
	    */
	    if(j != lb-1 && hr[j] == lb-1){
	      p = x + y*ncols + nrows*ncols*(j);
	      lg_label[p] = 0;
	      Seeds[k] = p;
	      k++;
	    }
	  }
	}
	else{ // it is background seed, then transfer for all layers
	  for(j = 0; j < nlayers ; j++){
	    p = x + y*ncols + nrows*ncols*(j);
	    lg_label[p] = 0;
	    Seeds[k] = p;
	    k++;
	  }
	}
      }
      Seeds[0] = k-1;

      /* ------------- Create LAYERED GRAPH -------------*/
      lg = gft::LayeredGraph::Create(nlayers, ncols*nrows);
      
      /*Create a graph image with 8-neighborhood*/
      sg = graph;
      //gft::SparseGraph::SuppressZeroWeightedArcs(sg);
      Wmax = sg->Wmax;

      /*Set cost of arcs for each layer*/
      for(i = 0; i < nlayers; i++){
	if(typeLayer[i] == 0){ /* 0 = normal*/
	  if(polLayer[i] == 0)
	    gft::LayeredGraph::SetArcs(lg, sg, i);
	  else{
	    sg_oriented = GetPolarityGraph(sg, cimg, img, polLayer[i]);
	    if(sg_oriented->Wmax > Wmax) Wmax = sg_oriented->Wmax;
	    gft::LayeredGraph::SetArcs(lg, sg_oriented, i);
	    gft::ImageGraph::Destroy(&sg_oriented);
	  }
	}
	else if(typeLayer[i] == 1){ /* 1 = GSC */
	  gft::sImageGraph *sg_GSC;
	  gft::sImage32 *P_sum;
	  int *Slb = NULL;
	  if(polLayer[i] == 0)
	    sg_GSC = gft::ImageGraph::Clone(sg);
	  else{
	    sg_GSC = GetPolarityGraph(sg, cimg, img, polLayer[i]);
	    if(sg_GSC->Wmax > Wmax) Wmax = sg_GSC->Wmax;
	  }
	  Slb = GetAllInternalSeedsByLabel(S, label, i+1, hr, nlayers);
	  //gft::SparseGraph::SuppressZeroWeightedArcs(sg_GSC);
	  P_sum = gft::ift::SC_Pred_fsum(sg_GSC, Slb, 0.1);
	  gft::ImageGraph::Orient2DigraphOuter(sg_GSC, P_sum);
	  gft::LayeredGraph::SetArcs(lg, sg_GSC, i);  
	  gft::Image32::Destroy(&P_sum);
	  gft::ImageGraph::Destroy(&sg_GSC);
	  gft::FreeIntArray(&Slb);
	}
      }

      
      for(i = 0; i < nlayers ; i++){
	if(hr[i] != -1){ //son and father
	  gft::LayeredGraph::SetArcs(lg, i, hr[i], ncols, 0.0, radius);
	  gft::LayeredGraph::SetArcs(lg, hr[i], i, ncols, (float)(Wmax+1), radius);
	}
      }
      

      //Second exclusion
      for(i = 0; i < nlayers-1; i++){
	for(j = i+1; j < nlayers; j++){
	  if( hr[i] == hr[j] ){ //same father -> brothers -> exclusion
	    gft::LayeredGraph::SetArcs(lg, i, j, ncols, (float)(Wmax+1), radius);
	    gft::LayeredGraph::SetArcs(lg, j, i, ncols, 0.0, radius);
	  }
	}
      }
      

      /* ------------- Executa a Hierarchical Layered OIFT -------------*/ 
      gft::ift::HL_OIFT(lg, Wmax, Seeds, lg_label, hr);

      
      /* ------------- OUTPUT RESULT-------------*/
      
      gft::sQueue *FIFO;
      FIFO = gft::Queue::Create(nlayers);
      int *depth = gft::AllocIntArray(nlayers);
      int dmax = 0;
      for(i = 0; i < nlayers ; i++){
	if(hr[i] == -1){
	  depth[i] = 0;
	  gft::Queue::Push(FIFO, i);
	}
	else
	  depth[i] = -1;
      }
      while(!gft::Queue::IsEmpty(FIFO)){
	j = gft::Queue::Pop(FIFO);
	for(i = 0; i < nlayers ; i++){
	  if(hr[i] == j){
	    depth[i] = depth[j] + 1;
	    if(depth[i] > dmax) dmax = depth[i];
	    gft::Queue::Push(FIFO, i);
	  }
	}
      }
      
      int sizeImg = ncols*nrows;
      
      gft::Image32::Set(label, 0);
      int d;
      for(d = 0; d <= dmax; d++){
	for(i = 0; i < nlayers; i++){
	  if(depth[i] == d){
	    for(j = 0; j < sizeImg; j++){
	      if(lg_label[j+ sizeImg*i] != 0)
		label->data[j] = i+1;
	    }
	  }
	}
      }

      gft::LayeredGraph::Destroy(&lg);
      
      gft::FreeIntArray(&typeLayer);
      gft::FreeIntArray(&polLayer);
      gft::FreeIntArray(&hr);
      gft::FreeIntArray(&lg_label);
      gft::FreeIntArray(&Seeds);
      gft::FreeIntArray(&depth);
      gft::Queue::Destroy(&FIFO);
    }


    //----------------------------------------

    
    void HL_OIFT(sLayeredGraph *lg,
		 int Wmax, int *S, int *L, int *hr){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,layer_p,layer_q,size_layer,tmp2=0;
      int w_pq=0, w_qp=0, tmp=0;
      int* value;
      sGraphNode *A;
      int exclusionCase = 0;
      int *i_inv = NULL;
      int i_inv_size;
      bool flag;

      A = lg->graph->nodes;
      p = 0;
      i_inv_size = A[p].outdegree;
      i_inv = gft::AllocIntArray(i_inv_size);
      flag = true;
      while(flag){
	flag = false;
	for(i = 0; i < A[p].outdegree && i < i_inv_size; i++){
	  q = A[p].adjList[i];
	  if(p == q){ flag = true; break; }
	  for(j = 0; j < A[q].outdegree; j++)
	    if(A[q].adjList[j] == p)
	      i_inv[i] = j;
	}
	p++;
      }
      
      size_layer = lg->nnodesperlayer;
      n = size_layer*lg->nlayers;
      
      //printf("Size Layer: %d\n",size_layer);
      //printf("N pixels in the Graph: %d\n",n);
      
      /*Initialization*/
      value = gft::AllocIntArray(n);
      //Q = gft::PQueue32::Create(Wmax*lg->nlayers+2,n,value);
      Q = gft::PQueue32::Create(Wmax+3,n,value);
      
      /*Insert in value*/
      for(p=0; p<n; p++){
	if(L[p]==NIL) value[p] = INT_MAX; //(int)floor(FLT_MAX+0.5);
	else          value[p] = 0;
      }
      
      /*Insert in PQueue*/
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  gft::PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(L[p]!=NIL)
	    gft::PQueue32::FastInsertElem(Q, p);	    
      }
      
      /*
	clock_t t;
	t = clock();
      */
      
      /*Starting OIFT*/
      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::FastRemoveMinFIFO(Q);
	
	for(i=0; i<A[p].outdegree; i++){
	  q = A[p].adjList[i];
	  //if(q != NIL){
	  if(Q->L.elem[q].color != BLACK){
	    //w_pq = ROUND(gft::Graph::GetArcWeight(lg->graph,p,q));
	    w_pq = ROUND(A[p].Warcs[i]);
	    //if(w_pq < 0) continue;
	    
	    layer_p = p/size_layer;
	    layer_q = q/size_layer;
	    exclusionCase = 0;
	    /* Analize each relation type*/
	    if (layer_p == layer_q){  /*SAME LAYER*/
	      // Get w_qp
	      //w_qp = ROUND(gft::Graph::GetArcWeight(lg->graph,q,p));
	      w_qp = ROUND(A[q].Warcs[i_inv[i]]);
	      exclusionCase = 0;
	      if(L[p] != 0)
		tmp = w_pq;
	      else 
		tmp = w_qp; 
	    }
	    /*INTER LAYERS: relations where defined*/
	    else if(layer_p != layer_q){  
	      // Get w_qp
	      if(w_pq == 0)
		w_qp = Wmax+1; //INT_MAX;
	      else
		w_qp = 0;
	      // inclusion
	      /*layer_p is "INCLUDED" in layer_q /OR/ layer_q is "INCLUDED" in layer_p */
	      if((hr[layer_p] == layer_q) || (hr[layer_q] == layer_p)){  
		exclusionCase = 0;
		if(L[p] != 0)
		  tmp = w_pq;
		else 
		  tmp = w_qp; 
	      }
	      // exclusion
	      /*layer_q is EXCLUDED from layer_p */
	      else if(hr[layer_p] == hr[layer_q]){  
		exclusionCase = 1;
		/* It was defined to change "inputs"/"outputs", in/from the higher layer for exclusion case*/
		if (layer_p < layer_q){ 
		  if(L[p] != 0)
		    tmp = w_qp;
		  else 
		    tmp = w_pq;
		}
		else{
		  if(L[p] != 0)
		    tmp = w_pq;
		  else 
		    tmp = w_qp;
		}
	      }
	    }
	    
	    if(tmp < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      
	      value[q] = tmp;
	      
	      if(exclusionCase == 1){ /*Treat exclusion case*/
		if(L[p] != 0)
		  L[q] = 0;
	      }
	      else{ /*exclusion == 0*/
		L[q] = L[p];
	      }
	      
	      gft::PQueue32::FastInsertElem(Q, q);
	    } 
	  }
	  //}
	}
      }

      /*
	t = clock() - t;
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("HLOIFT took %d seconds(t) to execute \n", ((int)t));
	printf("HLOIFT took %f seconds(t/C_P_S) to execute \n", time_taken);      
      */
      
      gft::PQueue32::Destroy(&Q);
      gft::FreeIntArray(&value);
      gft::FreeIntArray(&i_inv);
    }


    //--------------------------------------------------------
 
    void HL_OIFT_2(sLayeredGraph *lg,
		   int Wmax, int *S, int *L, int *hr){
      sPQueue32 *Q=NULL;
      int i,j,p,q,n,layer_p,layer_q,size_layer,tmp2=0;
      int w_pq=0, w_qp=0, tmp=0;
      int* value;
      sGraphNode *A;
      int exclusionCase = 0;

      A = lg->graph->nodes;
      size_layer = lg->nnodesperlayer;
      n = size_layer*lg->nlayers;

      /*Initialization*/
      value = gft::AllocIntArray(n);
      Q = gft::PQueue32::Create(Wmax+3,n,value);
      
      /*Insert in value*/
      for(p=0; p<n; p++){
	if(L[p]==NIL) value[p] = INT_MAX; //(int)floor(FLT_MAX+0.5);
	else          value[p] = 0;
      }
      
      /*Insert in PQueue*/
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  gft::PQueue32::FastInsertElem(Q, S[i]);
	}
      }
      else{
	for(p=0; p<n; p++)
	  if(L[p]!=NIL)
	    gft::PQueue32::FastInsertElem(Q, p);	    
      }
      
      /*Starting OIFT*/
      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::FastRemoveMinFIFO(Q);
	
	for(i=0; i<A[p].outdegree; i++){
	  q = A[p].adjList[i];
	  //if(q != NIL){
	  if(Q->L.elem[q].color != BLACK){
	    //w_pq = ROUND(gft::Graph::GetArcWeight(lg->graph,p,q));
	    w_pq = ROUND(A[p].Warcs[i]);
	    //if(w_pq < 0) continue;
	    
	    layer_p = p/size_layer;
	    layer_q = q/size_layer;
	    exclusionCase = 0;
	    /* Analize each relation type*/
	    if (layer_p == layer_q){  /*SAME LAYER*/
	      // Get w_qp
	      w_qp = ROUND(gft::Graph::GetArcWeight(lg->graph,q,p));
	      exclusionCase = 0;
	      if(L[p] != 0)
		tmp = w_pq;
	      else 
		tmp = w_qp; 
	    }
	    /*INTER LAYERS: relations where defined*/
	    else if(layer_p != layer_q){  
	      // Get w_qp
	      if(w_pq == 0)
		w_qp = Wmax+1; //INT_MAX;
	      else
		w_qp = 0;
	      // inclusion
	      /*layer_p is "INCLUDED" in layer_q /OR/ layer_q is "INCLUDED" in layer_p */
	      if((hr[layer_p] == layer_q) || (hr[layer_q] == layer_p)){  
		exclusionCase = 0;
		if(L[p] != 0)
		  tmp = w_pq;
		else 
		  tmp = w_qp; 
	      }
	      // exclusion
	      /*layer_q is EXCLUDED from layer_p */
	      else if(hr[layer_p] == hr[layer_q]){  
		exclusionCase = 1;
		/* It was defined to change "inputs"/"outputs", in/from the higher layer for exclusion case*/
		if (layer_p < layer_q){ 
		  if(L[p] != 0)
		    tmp = w_qp;
		  else 
		    tmp = w_pq;
		}
		else{
		  if(L[p] != 0)
		    tmp = w_pq;
		  else 
		    tmp = w_qp;
		}
	      }
	    }
	    
	    if(tmp < value[q]){
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::FastRemoveElem(Q, q);
	      
	      value[q] = tmp;
	      
	      if(exclusionCase == 1){ /*Treat exclusion case*/
		if(L[p] != 0)
		  L[q] = 0;
	      }
	      else{ /*exclusion == 0*/
		L[q] = L[p];
	      }
	      
	      gft::PQueue32::FastInsertElem(Q, q);
	    } 
	  }
	}
      }

      gft::PQueue32::Destroy(&Q);
      gft::FreeIntArray(&value);
    }
   


    //Compute watershed by fpeak from markers
    void IFT_fpeak(sImage32 *grad,
		   sAdjRel *A,
		   sImage32 *label){
      sPQueue32 *Q;
      sImage32 *cost; //*Rmin;
      int tmp, w, Wmax;
      int n,p,q,i,px,py,qx,qy;
      n = grad->n;
      Wmax = gft::Image32::GetMaxVal(grad);
      //Rmin = gft::Image32::RegMin(grad, A);
      cost = gft::Image32::Create(grad);
      Q = gft::PQueue32::Create(Wmax+2, n, cost->data);

      for(p = 0; p < n; p++){
	//if(Rmin->data[p] == 0){
	if(label->data[p] == NIL){
	  cost->data[p] = INT_MAX;
	  //label->data[p] = NIL;
	}
	else{
	  cost->data[p] = grad->data[p];
	  //label->data[p] = Rmin->data[p]-1;
	  PQueue32::FastInsertElem(Q, p);
	}
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	px = p%grad->ncols;
	py = p/grad->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(grad, qx,qy)){
	  if(qx >= 0 && qx < grad->ncols &&
	     qy >= 0 && qy < grad->nrows){
	    q = qx + qy*grad->ncols;
	      
	    if(Q->L.elem[q].color != BLACK){
	      w = grad->data[q];
	      tmp = MAX(cost->data[p], w);
	      
	      if(tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		cost->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      gft::Image32::Destroy(&cost);
      //gft::Image32::Destroy(&Rmin);
    }

    
    //Compute watershed by fwv from markers
    void IFT_fwv(sImage32 *grad,
		 sAdjRel *A,
		 sImage32 *label){
      sPQueue32 *Q;
      sImage32 *cost; //*Rmin;
      int tmp, Wmax;
      int n,p,q,i,px,py,qx,qy;
      n = grad->n;
      Wmax = gft::Image32::GetMaxVal(grad);
      //Rmin = gft::Image32::RegMin(grad, A);
      cost = gft::Image32::Create(grad);
      Q = gft::PQueue32::Create(Wmax+2, n, cost->data);

      for(p = 0; p < n; p++){
	//if(Rmin->data[p] == 0){
	if(label->data[p] == NIL){
	  cost->data[p] = INT_MAX;
	  //label->data[p] = NIL;
	}
	else{
	  cost->data[p] = grad->data[p];
	  //label->data[p] = Rmin->data[p]-1;
	  PQueue32::FastInsertElem(Q, p);
	}
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	px = p%grad->ncols;
	py = p/grad->ncols;
	for(i=1; i < A->n; i++){
	  qx = px + A->dx[i];
	  qy = py + A->dy[i];
	  //if(gft::Image32::IsValidPixel(grad, qx,qy)){
	  if(qx >= 0 && qx < grad->ncols &&
	     qy >= 0 && qy < grad->nrows){
	      q = qx + qy*grad->ncols;
	      
	    if(Q->L.elem[q].color != BLACK){
	      tmp = grad->data[q];
	      
	      if(tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		cost->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      gft::Image32::Destroy(&cost);
      //gft::Image32::Destroy(&Rmin);
    }



    //-------------------------------
    

  } /*end ift namespace*/
} /*end gft namespace*/


