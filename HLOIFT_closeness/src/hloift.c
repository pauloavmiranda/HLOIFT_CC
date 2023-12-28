
#include "hloift.h"


gft::sLayeredGraph *HL_OIFT_CreateGraph(gft::sImageGraph *graph,
					gft::sCImage *cimg,
					gft::sImage32 *img,
					float radius,
					char *hierarchy,
					int **p_hr,
					int *p_Wmax,
					gft::sImage32 **p_P,
					int x, int y){
  gft::sLayeredGraph *lg;
  gft::sImageGraph *sg, *sg_oriented;
  FILE *fpHierarchy;
  int nlayers, typeL, nRelations, relationT, layer1, layer2;
  int i, j, n, Wmax, pol;
  int ncols, nrows;
  //---------------------
  gft::sImageGraph *sg_GSC;
  gft::sImage32 *P = NULL;
  bool computeGSC = false;
  int Sobj[2];
  
  ncols = graph->ncols;
  nrows = graph->nrows;
  gft::ImageGraph::ChangeType(graph, DISSIMILARITY);

  fpHierarchy = fopen((char *)hierarchy, "r");
  if (!fpHierarchy){
    gft::Error((char *)MSG1,(char *)"Error opening file of Hierarchy Relations!");
  }

  /*------------- Set number of layers -------------*/
  fscanf(fpHierarchy," %d\n",&nlayers);
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
    //printf("typeL: %d, pol: %d\n", typeL, pol);
    if(typeL == 1)
      computeGSC = true;
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
  
  /* ------------- Create LAYERED GRAPH -------------*/
  lg = gft::LayeredGraph::Create(nlayers, ncols*nrows);

  if(computeGSC){
    Sobj[0] = 1;
    Sobj[1] = x + y*graph->ncols;
    P = gft::ift::SC_Pred_fsum(graph, Sobj, 0.1);
    *p_P = P;
  }
  
  /*Create a graph image with 8-neighborhood*/
  sg = graph;
  Wmax = sg->Wmax;
  
  /*Set cost of arcs for each layer*/
  for(i = 0; i < nlayers; i++){
    if(typeLayer[i] == 0){ /* 0 = normal*/
      if(polLayer[i] == 0)
	gft::LayeredGraph::SetArcs(lg, sg, i);
      else{
	sg_oriented = gft::ift::GetPolarityGraph(sg, cimg, img, polLayer[i]);
	if(sg_oriented->Wmax > Wmax) Wmax = sg_oriented->Wmax;
	gft::LayeredGraph::SetArcs(lg, sg_oriented, i);
	gft::ImageGraph::Destroy(&sg_oriented);
      }
    }
    else if(typeLayer[i] == 1){ // 1 = GSC
      if(polLayer[i] == 0)
	sg_GSC = gft::ImageGraph::Clone(sg);
      else{
	sg_GSC = gft::ift::GetPolarityGraph(sg, cimg, img, polLayer[i]);
	if(sg_GSC->Wmax > Wmax) Wmax = sg_GSC->Wmax;
      }
      if(P != NULL)
	gft::ImageGraph::Orient2DigraphOuter(sg_GSC, P);
      gft::LayeredGraph::SetArcs(lg, sg_GSC, i);  
      gft::ImageGraph::Destroy(&sg_GSC);
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
  
  gft::FreeIntArray(&typeLayer);
  gft::FreeIntArray(&polLayer);
  *p_hr = hr;
  *p_Wmax = Wmax;
  return lg;
}

//----------------------------------------

int *HL_OIFT_CreateSeeds(gft::sLayeredGraph *lg,
			 int *hr,
			 int *S,
			 gft::sImage32 *label,
			 int **p_lg_label){
  int *lg_label;
  int *Seeds;
  int n,i,j,k,x,y,s,p,lb,nSeeds;
  
  n = label->n*lg->nlayers;
  lg_label = gft::AllocIntArray(n);
  for(i = 0; i < n; i++){
    lg_label[i] = NIL;
  }
  
  nSeeds = S[0];
  Seeds = gft::AllocIntArray(nSeeds*lg->nlayers + 1);
  
  k = 1;
  for(i = 1; i <= S[0]; i++){
    s = S[i];
    x = s%label->ncols;
    y = s/label->ncols;
    lb = label->data[s];
    
    if(lb != 0){
      p = x + y*label->ncols + label->n*(lb-1); 
      lg_label[p] = 1;
      Seeds[k] = p;
      k++;
      /*Essa parte faltou no algoritmo do artigo:*/
      for(j = 0; j < lg->nlayers; j++){
	if(j != lb-1 && hr[j] == lb-1){
	  p = x + y*label->ncols + label->n*(j);
	  lg_label[p] = 0;
	  Seeds[k] = p;
	  k++;
	}
      }
    }
    else{ // it is background seed, then transfer for all layers
      for(j = 0; j < lg->nlayers ; j++){
	p = x + y*label->ncols + label->n*(j);
	lg_label[p] = 0;
	Seeds[k] = p;
	k++;
      }
    }
  }
  Seeds[0] = k-1;
  *p_lg_label = lg_label;
  return Seeds;
}


//----------------------------------------

void HL_OIFT_Segmentation(gft::sLayeredGraph *lg,
			  int Wmax,
			  int *hr,
			  int *S,
			  gft::sImage32 *label){
  gft::sImage32 *tmp = NULL;
  int *lg_label = NULL;
  int *Seeds = NULL;
  int i,j;

  Seeds = HL_OIFT_CreateSeeds(lg, hr, S, label, &lg_label);

  // ------------- Executa a Hierarchical Layered OIFT -------------
  HL_OIFT_Algorithm(lg, Wmax, Seeds, lg_label, hr);
  
  
  // ------------- OUTPUT RESULT-------------
      
  gft::sQueue *FIFO;
  FIFO = gft::Queue::Create(lg->nlayers);
  int *depth = gft::AllocIntArray(lg->nlayers);
  int dmax = 0;
  for(i = 0; i < lg->nlayers; i++){
    if(hr[i] == -1){
      depth[i] = 0;
      gft::Queue::Push(FIFO, i);
    }
    else
      depth[i] = -1;
  }
  while(!gft::Queue::IsEmpty(FIFO)){
    j = gft::Queue::Pop(FIFO);
    for(i = 0; i < lg->nlayers; i++){
      if(hr[i] == j){
	depth[i] = depth[j] + 1;
	if(depth[i] > dmax) dmax = depth[i];
	gft::Queue::Push(FIFO, i);
      }
    }
  }
  
  int sizeImg = label->n;
  
  gft::Image32::Set(label, 0);
  int d;
  for(d = 0; d <= dmax; d++){
    for(i = 0; i < lg->nlayers; i++){
      if(depth[i] == d){
	for(j = 0; j < sizeImg; j++){
	  if(lg_label[j+ sizeImg*i] != 0)
	    label->data[j] = i+1;
	}
      }
    }
  }
  
  gft::FreeIntArray(&depth);
  gft::Queue::Destroy(&FIFO);
  gft::FreeIntArray(&Seeds);
  gft::FreeIntArray(&lg_label);
}



void HL_OIFT_Algorithm(gft::sLayeredGraph *lg,
		       int Wmax, int *S, int *L, int *hr){
  gft::sPQueue32 *Q=NULL;
  int i,p,q,n,layer_p,layer_q,size_layer;
  int w_pq=0, w_qp=0, tmp=0;
  int *value;
  gft::sGraphNode *A;
  int exclusionCase = 0;
  int *i_inv = NULL;

  A = lg->graph->nodes;
  i_inv = IndexesOfAntiparallelArcs(lg);
      
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




//----------------------------------------
gft::sImage32 *Get_Geodesic_Distance(gft::sImage32 *pred,
				     gft::sAdjRel *A){
  gft::sQueue *Qfifo=NULL;
  gft::sImage32 *cost;
  gft::Pixel u,v;
  int n,p,q,i;
  int *Dpq;
  
  Dpq = gft::AllocIntArray(A->n);
  for(i = 1; i < A->n; i++){
    Dpq[i] = ROUND(10.0*sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]));
  }
  
  n = pred->n;
  cost  = gft::Image32::Create(pred->ncols, pred->nrows);
  gft::Image32::Set(cost, INT_MAX);
  Qfifo = gft::Queue::Create(n);
      
  for(p = 0; p < n; p++){
    if(pred->data[p] == NIL){
      cost->data[p] = 0;
      gft::Queue::Push(Qfifo, p);
    }
  }
  
  while(!gft::Queue::IsEmpty(Qfifo)){
    p = gft::Queue::Pop(Qfifo);
    u.x = p%pred->ncols; 
    u.y = p/pred->ncols; 
    
    for(i = 1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if(v.x >= 0 && v.x < pred->ncols &&
	 v.y >= 0 && v.y < pred->nrows){
	q = v.x + pred->ncols*v.y;
	if(p == pred->data[q]){
	  cost->data[q] = cost->data[p] + Dpq[i];
	  gft::Queue::Push(Qfifo, q);
	}
      }
    }
  }
  gft::FreeIntArray(&Dpq);
  gft::Queue::Destroy(&Qfifo);
  return cost;
}




void HL_OIFT_Closeness_Segmentation(gft::sLayeredGraph *lg,
				    int Wmax,
				    int *hr,
				    int *S,
				    int cr_L,
				    gft::sImage32 *P,
				    gft::sImage32 *label,
				    gft::sAdjRel *A){
  gft::sImage32 *C = NULL;
  int *lg_label = NULL;
  int *Seeds = NULL;
  int i,j,l;

  if(lg->nlayers > 2){
    printf("We currently only support two layers.\n");
    exit(1);
  }

  Seeds = HL_OIFT_CreateSeeds(lg, hr, S, label, &lg_label);

  l = lg->nlayers - 1;

  C = Get_Geodesic_Distance(P, A);
  cr_L *= 10;
  // ------------- Executa a Hierarchical Layered OIFT -----------
  HL_OIFT_Closeness_Algorithm(lg, Wmax, Seeds, lg_label,
			      hr, cr_L, C, P, A);
  
  // ------------- OUTPUT RESULT-------------
      
  gft::sQueue *FIFO;
  FIFO = gft::Queue::Create(lg->nlayers);
  int *depth = gft::AllocIntArray(lg->nlayers);
  int dmax = 0;
  for(i = 0; i < lg->nlayers; i++){
    if(hr[i] == -1){
      depth[i] = 0;
      gft::Queue::Push(FIFO, i);
    }
    else
      depth[i] = -1;
  }
  while(!gft::Queue::IsEmpty(FIFO)){
    j = gft::Queue::Pop(FIFO);
    for(i = 0; i < lg->nlayers; i++){
      if(hr[i] == j){
	depth[i] = depth[j] + 1;
	if(depth[i] > dmax) dmax = depth[i];
	gft::Queue::Push(FIFO, i);
      }
    }
  }
  
  int sizeImg = label->n;
  
  gft::Image32::Set(label, 0);
  int d;
  for(d = 0; d <= dmax; d++){
    for(i = 0; i < lg->nlayers; i++){
      if(depth[i] == d){
	for(j = 0; j < sizeImg; j++){
	  if(lg_label[j+ sizeImg*i] != 0)
	    label->data[j] = i+1;
	}
      }
    }
  }

  gft::Image32::Destroy(&C);
  gft::FreeIntArray(&depth);
  gft::Queue::Destroy(&FIFO);
  gft::FreeIntArray(&Seeds);
  gft::FreeIntArray(&lg_label);
}


//----------------------------------------


int *IndexesOfAntiparallelArcs(gft::sLayeredGraph *lg){
  gft::sGraphNode *A;
  int *i_inv = NULL;
  int i_inv_size, p,q,i,j, layer_p,layer_q;
  bool flag;

  A = lg->graph->nodes;
  p = 0;
  i_inv_size = 0;
  for(i = 0; i < A[p].outdegree; i++){
    q = A[p].adjList[i];
    layer_p = p/lg->nnodesperlayer;
    layer_q = q/lg->nnodesperlayer;
    if(layer_p != layer_q)
      break;
    i_inv_size++;
  }

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
  return i_inv;
}

//----------------------------------------

void Prune_Pixels_For_Closeness(int r,
				int cr_L,
				gft::sImage32 *C,
				gft::sImage32 *P,
				int *L,
				gft::sPQueue32 *Q,
				gft::sAdjRel *A,
				int *T){
  gft::Pixel u,v;
  int i, p, q, qq, n, size_layer, Cpr;
  size_layer = P->n;
  Cpr = C->data[P->data[r]];
  n = 1;
  T[0] = r;
  while(n > 0){
    n--;
    p = T[n];
    u.x = p % P->ncols;
    u.y = p / P->ncols;    
    for(i = 1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if(v.x >= 0 && v.x < P->ncols &&
	 v.y >= 0 && v.y < P->nrows){
	q = v.x + v.y*P->ncols;
	if(p == P->data[q]){
	  qq = q + 0*size_layer;
	  if(C->data[q] - Cpr > cr_L &&
	     Q->L.elem[qq].color != BLACK){
	    L[qq] = 0;
	    if(Q->L.elem[qq].color == GRAY)
	      gft::PQueue32::FastRemoveElem(Q, qq);
	    Q->L.value[qq] = 0;
	    gft::PQueue32::FastInsertElem(Q, qq);	    
	  }
	  else if(!(Q->L.elem[qq].color == BLACK &&
		    L[qq] == 0)){ //&& C->data[q] - Cpr <= cr_L){
	    T[n] = q;
	    n++;
	  }
	}
      }
    }
  }
}



void Conquest_Support_Pixel_For_Closeness(int v,
					  int cr_L,
					  gft::sImage32 *C,
					  gft::sImage32 *P,
					  int *L,
					  gft::sPQueue32 *Q){
  int t, p, tt, size_layer;
  size_layer = P->n;
  t = v;
  p = P->data[v];
  while(!(Q->L.elem[p + size_layer].color == BLACK &&
	  L[p + size_layer] == 1)  &&
	C->data[v]-C->data[p] <= cr_L){
    t = p;
    p = P->data[p];
  }
  if(C->data[v]-C->data[p] > cr_L){
    tt = t + size_layer;
    if(Q->L.elem[tt].color != BLACK){
      L[tt] = 1;
      if(Q->L.elem[tt].color == GRAY)
	gft::PQueue32::FastRemoveElem(Q, tt);
      Q->L.value[tt] = 0;
      gft::PQueue32::FastInsertElem(Q, tt);
    }
  }
}


//----------------------------------------


void HL_OIFT_Closeness_Algorithm(gft::sLayeredGraph *lg,
				 int Wmax, int *S, int *L, int *hr,
				 int cr_L,
				 gft::sImage32 *C,
				 gft::sImage32 *P,
				 gft::sAdjRel *A){
  gft::sPQueue32 *Q=NULL;
  int i,p,pp,q,n,layer_p,layer_q,size_layer;
  int w_pq=0, w_qp=0, tmp=0;
  int *value;
  gft::sGraphNode *GN;
  int exclusionCase = 0;
  int *i_inv = NULL;
  int t,nsteps = 0;
  int *Stack;
  //------------------------
  //int total_prune = 0;
  
  Stack = (int *)malloc(sizeof(int)*P->n);
  GN = lg->graph->nodes;
  i_inv = IndexesOfAntiparallelArcs(lg);
      
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

    layer_p = p/size_layer;
    //------------------------------------
    if(L[p] > 0){
      if(layer_p == 0){
	Conquest_Support_Pixel_For_Closeness(p%size_layer,
					     cr_L, C, P, L, Q);
      }
    }
    else{ //L[p] == 0
      if(layer_p == 1){
	pp = p%size_layer;
	t = P->data[pp] + 1*size_layer;
	if(!(L[t] == 0 && Q->L.elem[t].color == BLACK)){
	  //total_prune += 1;
	  Prune_Pixels_For_Closeness(pp, cr_L, C, P,
				     L, Q, A, Stack);
	}
      }
    }
    //------------------------------------
    
    for(i=0; i<GN[p].outdegree; i++){
      q = GN[p].adjList[i];
      //if(q != NIL){
      if(Q->L.elem[q].color != BLACK){
	//w_pq = ROUND(gft::Graph::GetArcWeight(lg->graph,p,q));
	w_pq = ROUND(GN[p].Warcs[i]);
	//if(w_pq < 0) continue;
	
	layer_q = q/size_layer;
	exclusionCase = 0;
	/* Analize each relation type*/
	if (layer_p == layer_q){  /*SAME LAYER*/
	  // Get w_qp
	  //w_qp = ROUND(gft::Graph::GetArcWeight(lg->graph,q,p));
	  w_qp = ROUND(GN[q].Warcs[i_inv[i]]);
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
  free(Stack);
  gft::PQueue32::Destroy(&Q);
  gft::FreeIntArray(&value);
  gft::FreeIntArray(&i_inv);

  //printf("total_prune: %d\n", total_prune);
}




