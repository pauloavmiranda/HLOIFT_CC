#include "gft_heap32fiif_lex.h"

namespace gft{
  namespace Heap32fiif_lex{


    void GoUp_MaxPolicy(sHeap32fiif_lex *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      float c1 = H->cost1[p];
      int   c2 = H->cost2[p];
      int   c3 = H->cost3[p];
      float c4 = H->cost4[p];
      while((i > 1)&&
	    ((H->cost1[H->pixel[j]] < c1)||
	     (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] < c2)||
	     (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] < c3)||
	     (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] == c3 && H->cost4[H->pixel[j]] < c4)
	     )){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }


    void GoUp_MinPolicy(sHeap32fiif_lex *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      float c1 = H->cost1[p];
      int   c2 = H->cost2[p];
      int   c3 = H->cost3[p];
      float c4 = H->cost4[p];
      while((i > 1)&&
	    ((H->cost1[H->pixel[j]] > c1)||
	     (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] > c2)||
	     (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] > c3)||
	     (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] == c3 && H->cost4[H->pixel[j]] > c4)
	     )){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }    

    

    void GoDown_MaxPolicy(sHeap32fiif_lex *H, int i) {
      int j, left, right, p = H->pixel[i];
      float c1 = H->cost1[p];
      int   c2 = H->cost2[p];
      int   c3 = H->cost3[p];
      float c4 = H->cost4[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(right <= H->last){
	  if((H->cost1[H->pixel[right]] > H->cost1[H->pixel[left]])||
	     (H->cost1[H->pixel[right]] == H->cost1[H->pixel[left]] && H->cost2[H->pixel[right]] > H->cost2[H->pixel[left]])||
	     (H->cost1[H->pixel[right]] == H->cost1[H->pixel[left]] && H->cost2[H->pixel[right]] == H->cost2[H->pixel[left]] && H->cost3[H->pixel[right]] > H->cost3[H->pixel[left]])||
	     (H->cost1[H->pixel[right]] == H->cost1[H->pixel[left]] && H->cost2[H->pixel[right]] == H->cost2[H->pixel[left]] && H->cost3[H->pixel[right]] == H->cost3[H->pixel[left]] && H->cost4[H->pixel[right]] > H->cost4[H->pixel[left]])
	     )
	    j = right;
	  else
	    j = left;
	}
	else if(left <= H->last)
	  j = left;
	else
	  break;
	
	if((H->cost1[H->pixel[j]] > c1)||
	   (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] > c2)||
	   (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] > c3)||
	   (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] == c3 && H->cost4[H->pixel[j]] > c4)
	   ){
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }
    

    void GoDown_MinPolicy(sHeap32fiif_lex *H, int i) {
      int j, left, right, p = H->pixel[i];
      float c1 = H->cost1[p];
      int   c2 = H->cost2[p];
      int   c3 = H->cost3[p];
      float c4 = H->cost4[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(right <= H->last){
	  if((H->cost1[H->pixel[right]] < H->cost1[H->pixel[left]])||
	     (H->cost1[H->pixel[right]] == H->cost1[H->pixel[left]] && H->cost2[H->pixel[right]] < H->cost2[H->pixel[left]])||
	     (H->cost1[H->pixel[right]] == H->cost1[H->pixel[left]] && H->cost2[H->pixel[right]] == H->cost2[H->pixel[left]] && H->cost3[H->pixel[right]] < H->cost3[H->pixel[left]])||
	     (H->cost1[H->pixel[right]] == H->cost1[H->pixel[left]] && H->cost2[H->pixel[right]] == H->cost2[H->pixel[left]] && H->cost3[H->pixel[right]] == H->cost3[H->pixel[left]] && H->cost4[H->pixel[right]] < H->cost4[H->pixel[left]])
	     )
	    j = right;
	  else
	    j = left;
	}
	else if(left <= H->last)
	  j = left;
	else
	  break;
	
	if((H->cost1[H->pixel[j]] < c1)||
	   (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] < c2)||
	   (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] < c3)||
	   (H->cost1[H->pixel[j]] == c1 && H->cost2[H->pixel[j]] == c2 && H->cost3[H->pixel[j]] == c3 && H->cost4[H->pixel[j]] < c4)
	   ){
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }

    
    
    char IsFull(sHeap32fiif_lex *H) {
      if (H->last == H->n)
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(sHeap32fiif_lex *H) {
      if (H->last == 0)
	return 1;
      else
	return 0;
    }

    
    sHeap32fiif_lex *Create(int n,
			    float *cost1, int   *cost2,
			    int   *cost3, float *cost4) {
      sHeap32fiif_lex *H = NULL;
      int i;
      
      if (cost1 == NULL || cost2 == NULL || cost3 == NULL || cost4 == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap32fiif_lex::Create");
	return NULL;
      }
      
      H = (sHeap32fiif_lex *) malloc(sizeof(sHeap32fiif_lex));
      if (H != NULL) {
	H->n       = n;
	H->cost1   = cost1;
	H->cost2   = cost2;
	H->cost3   = cost3;
	H->cost4   = cost4;
	H->color   = (char *) malloc(sizeof(char) * n);
	H->pixel   = (int *) malloc(sizeof(int) * (n+1));
	H->pos     = (int *) malloc(sizeof(int) * n);
	H->last    = 0;
	if (H->color == NULL || H->pos == NULL || H->pixel == NULL)
	  gft::Error((char *)MSG1,(char *)"Heap32fiif_lex::Create");
	for (i = 0; i < H->n; i++) {
	  H->color[i] = WHITE;
	  //H->pos[i]   = -1;
	  //H->pixel[i] = -1;
	}
	//H->pixel[n] = -1;
      }
      else
	gft::Error((char *)MSG1,(char *)"Heap32fiif_lex::Create");

      return H;
    }

    
    void Destroy(sHeap32fiif_lex **H) {
      sHeap32fiif_lex *aux = *H;
      if (aux != NULL) {
	if (aux->pixel != NULL) free(aux->pixel);
	if (aux->color != NULL) free(aux->color);
	if (aux->pos != NULL)   free(aux->pos);
	free(aux);
	*H = NULL;
      }
    }

    
    void Insert_MaxPolicy(sHeap32fiif_lex *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MaxPolicy(H, H->last);
    }

    
    void Insert_MinPolicy(sHeap32fiif_lex *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MinPolicy(H, H->last);
    }

    
    void Remove_MaxPolicy(sHeap32fiif_lex *H, int *pixel) {
      *pixel = H->pixel[1];
      //H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	//H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else{ //else if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	//H->pixel[H->last] = -1;
	H->last--;
	GoDown_MaxPolicy(H, 1);
      }
    }


    void Remove_MinPolicy(sHeap32fiif_lex *H, int *pixel) {
      *pixel = H->pixel[1];
      //H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	//H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else{ //else if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	//H->pixel[H->last] = -1;
	H->last--;
	GoDown_MinPolicy(H, 1);
      }
    }


    void Get_MaxPolicy(sHeap32fiif_lex *H, int *pixel){
      *pixel = H->pixel[1];
    }

    
    void Get_MinPolicy(sHeap32fiif_lex *H, int *pixel){
      *pixel = H->pixel[1];
    }
    

    void Update_MaxPolicy(sHeap32fiif_lex *H, int p,
			  float value1, int   value2,
			  int   value3, float value4){
      H->cost1[p] = value1;
      H->cost2[p] = value2;
      H->cost3[p] = value3;
      H->cost4[p] = value4;
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MaxPolicy(H, p);
      else
	GoUp_MaxPolicy(H, H->pos[p]);
    }


    void Update_MinPolicy(sHeap32fiif_lex *H, int p,
			  float value1, int   value2,
			  int   value3, float value4){
      H->cost1[p] = value1;
      H->cost2[p] = value2;
      H->cost3[p] = value3;
      H->cost4[p] = value4;
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MinPolicy(H, p);
      else
	GoUp_MinPolicy(H, H->pos[p]);
    }
    
    
    void Reset(sHeap32fiif_lex *H){
      int i;
      
      for (i=0; i < H->n; i++) {
	H->color[i] = WHITE;
	//H->pos[i]   = -1;
	//H->pixel[i] = -1;
      }
      //H->pixel[H->n] = -1;
      H->last = 0;
    }



    void Delete_MaxPolicy(sHeap32fiif_lex *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      H->color[pixel] = WHITE;
      q = H->pixel[H->last];
      if(pixel != q){
	H->pixel[H->pos[pixel]] = q;
	H->pos[q] = H->pos[pixel];
	//H->pixel[H->last] = -1;
	H->last--;
	//H->pos[pixel] = -1;
	if((H->cost1[pixel] > H->cost1[q])||
	   (H->cost1[pixel] == H->cost1[q] && H->cost2[pixel] > H->cost2[q])||
	   (H->cost1[pixel] == H->cost1[q] && H->cost2[pixel] == H->cost2[q] && H->cost3[pixel] > H->cost3[q])||
	   (H->cost1[pixel] == H->cost1[q] && H->cost2[pixel] == H->cost2[q] && H->cost3[pixel] == H->cost3[q] && H->cost4[pixel] > H->cost4[q])
	   )
	  GoDown_MaxPolicy(H, H->pos[q]);
	else 
	  GoUp_MaxPolicy(H, H->pos[q]);
      }
      else{
	//H->pixel[H->last] = -1;
	H->last--;
	//H->pos[pixel] = -1;
      }
    }
    

    void Delete_MinPolicy(sHeap32fiif_lex *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      H->color[pixel] = WHITE;
      q = H->pixel[H->last];
      if(pixel != q){
	H->pixel[H->pos[pixel]] = q;
	H->pos[q] = H->pos[pixel];
	//H->pixel[H->last] = -1;
	H->last--;
	//H->pos[pixel] = -1;
	if((H->cost1[pixel] < H->cost1[q])||
	   (H->cost1[pixel] == H->cost1[q] && H->cost2[pixel] < H->cost2[q])||
	   (H->cost1[pixel] == H->cost1[q] && H->cost2[pixel] == H->cost2[q] && H->cost3[pixel] < H->cost3[q])||
	   (H->cost1[pixel] == H->cost1[q] && H->cost2[pixel] == H->cost2[q] && H->cost3[pixel] == H->cost3[q] && H->cost4[pixel] < H->cost4[q])
	   )
	  GoDown_MinPolicy(H, H->pos[q]);
	else 
	  GoUp_MinPolicy(H, H->pos[q]);
      }
      else{
	//H->pixel[H->last] = -1;
	H->last--;
	//H->pos[pixel] = -1;
      }
    }

    

  } /*end Heap32fiif_lex namespace*/
} /*end gft namespace*/


