
LIB=./lib
INCLUDE=./include
BIN=./bin
SRC=./src
OBJ=./obj

#Compiladores
CXX=g++

#-march=native
#FLAGS= -Wall -msse -fopenmp
#FLAGS= -Wall -O3 -msse -fopenmp
#FLAGS= -Wall -g
#FLAGS= -Wall -g -O3
FLAGS= -Wall -O3
#FLAGS= -Wall

LINKS= -lpthread -lz -lm 


all: libgft
	@echo
	@echo "libgft.a built..."
	@echo "DONE."

libgft: $(LIB)/libgft.a

GFTLIB  = -L$(LIB) -lgft
GFTFLAGS  = -I$(INCLUDE)


$(LIB)/libgft.a: \
$(OBJ)/gft_common.o \
$(OBJ)/gft_adjrel.o \
$(OBJ)/gft_image32.o \
$(OBJ)/gft_image32f.o \
$(OBJ)/gft_image64.o \
$(OBJ)/gft_cimage.o \
$(OBJ)/gft_cimage32f.o \
$(OBJ)/gft_morphology.o \
$(OBJ)/gft_queue.o \
$(OBJ)/gft_gpqueue_by_Falcao.o \
$(OBJ)/gft_graph.o \
$(OBJ)/gft_imagegraph.o \
$(OBJ)/gft_imagegraphpx.o \
$(OBJ)/gft_layeredgraph.o \
$(OBJ)/gft_heap.o \
$(OBJ)/gft_heap32.o \
$(OBJ)/gft_heap64f.o \
$(OBJ)/gft_heap_lex.o \
$(OBJ)/gft_heap32fi_lex.o \
$(OBJ)/gft_heap32fif_lex.o \
$(OBJ)/gft_heap32fifi_lex.o \
$(OBJ)/gft_heap32fiif_lex.o \
$(OBJ)/gft_heappair.o \
$(OBJ)/gft_ift.o \
$(OBJ)/gft_string.o \
$(OBJ)/gft_arraylist.o \
$(OBJ)/gft_filelist.o \
$(OBJ)/gft_color.o \
$(OBJ)/gft_highlight.o \
$(OBJ)/gft_evaluation.o \
$(OBJ)/gft_matrix.o \
$(OBJ)/gft_bmap.o \
$(OBJ)/gft_set.o \
$(OBJ)/gft_stack.o \
$(OBJ)/gft_analysis.o \
$(OBJ)/gft_pqueue32.o \
$(OBJ)/gft_pqueue16.o \
$(OBJ)/gft_hashtable.o \
$(OBJ)/gft_curve.o \
$(OBJ)/gft_radiometric.o \
$(OBJ)/gft_filtering.o \
$(OBJ)/gft_marker.o \
$(OBJ)/gft_attributelist.o \
$(OBJ)/gft_configfile.o \
$(OBJ)/gft_geometry.o \
$(OBJ)/gft_draw.o
	ar csr $(LIB)/libgft.a \
$(OBJ)/gft_common.o \
$(OBJ)/gft_adjrel.o \
$(OBJ)/gft_image32.o \
$(OBJ)/gft_image32f.o \
$(OBJ)/gft_image64.o \
$(OBJ)/gft_cimage.o \
$(OBJ)/gft_cimage32f.o \
$(OBJ)/gft_morphology.o \
$(OBJ)/gft_queue.o \
$(OBJ)/gft_gpqueue_by_Falcao.o \
$(OBJ)/gft_graph.o \
$(OBJ)/gft_imagegraph.o \
$(OBJ)/gft_imagegraphpx.o \
$(OBJ)/gft_layeredgraph.o \
$(OBJ)/gft_heap.o \
$(OBJ)/gft_heap32.o \
$(OBJ)/gft_heap64f.o \
$(OBJ)/gft_heap_lex.o \
$(OBJ)/gft_heap32fi_lex.o \
$(OBJ)/gft_heap32fif_lex.o \
$(OBJ)/gft_heap32fifi_lex.o \
$(OBJ)/gft_heap32fiif_lex.o \
$(OBJ)/gft_heappair.o \
$(OBJ)/gft_ift.o \
$(OBJ)/gft_string.o \
$(OBJ)/gft_arraylist.o \
$(OBJ)/gft_filelist.o \
$(OBJ)/gft_color.o \
$(OBJ)/gft_highlight.o \
$(OBJ)/gft_evaluation.o \
$(OBJ)/gft_matrix.o \
$(OBJ)/gft_bmap.o \
$(OBJ)/gft_set.o \
$(OBJ)/gft_stack.o \
$(OBJ)/gft_analysis.o \
$(OBJ)/gft_pqueue32.o \
$(OBJ)/gft_pqueue16.o \
$(OBJ)/gft_hashtable.o \
$(OBJ)/gft_curve.o \
$(OBJ)/gft_radiometric.o \
$(OBJ)/gft_filtering.o \
$(OBJ)/gft_marker.o \
$(OBJ)/gft_attributelist.o \
$(OBJ)/gft_configfile.o \
$(OBJ)/gft_geometry.o \
$(OBJ)/gft_draw.o


$(OBJ)/gft_common.o: $(SRC)/gft_common.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_common.cpp -I$(INCLUDE) -o $(OBJ)/gft_common.o

$(OBJ)/gft_adjrel.o: $(SRC)/gft_adjrel.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_adjrel.cpp -I$(INCLUDE) -o $(OBJ)/gft_adjrel.o

$(OBJ)/gft_image32.o: $(SRC)/gft_image32.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_image32.cpp -I$(INCLUDE) -o $(OBJ)/gft_image32.o

$(OBJ)/gft_image32f.o: $(SRC)/gft_image32f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_image32f.cpp -I$(INCLUDE) -o $(OBJ)/gft_image32f.o

$(OBJ)/gft_image64.o: $(SRC)/gft_image64.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_image64.cpp -I$(INCLUDE) -o $(OBJ)/gft_image64.o

$(OBJ)/gft_cimage.o: $(SRC)/gft_cimage.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_cimage.cpp -I$(INCLUDE) -o $(OBJ)/gft_cimage.o

$(OBJ)/gft_cimage32f.o: $(SRC)/gft_cimage32f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_cimage32f.cpp -I$(INCLUDE) -o $(OBJ)/gft_cimage32f.o

$(OBJ)/gft_morphology.o: $(SRC)/gft_morphology.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_morphology.cpp -I$(INCLUDE) -o $(OBJ)/gft_morphology.o

$(OBJ)/gft_queue.o: $(SRC)/gft_queue.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_queue.cpp -I$(INCLUDE) -o $(OBJ)/gft_queue.o

$(OBJ)/gft_gpqueue_by_Falcao.o: $(SRC)/gft_gpqueue_by_Falcao.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_gpqueue_by_Falcao.cpp -I$(INCLUDE) -o $(OBJ)/gft_gpqueue_by_Falcao.o

$(OBJ)/gft_graph.o: $(SRC)/gft_graph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_graph.cpp -I$(INCLUDE) -o $(OBJ)/gft_graph.o

$(OBJ)/gft_imagegraph.o: $(SRC)/gft_imagegraph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_imagegraph.cpp -I$(INCLUDE) -o $(OBJ)/gft_imagegraph.o

$(OBJ)/gft_imagegraphpx.o: $(SRC)/gft_imagegraphpx.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_imagegraphpx.cpp -I$(INCLUDE) -o $(OBJ)/gft_imagegraphpx.o

$(OBJ)/gft_layeredgraph.o: $(SRC)/gft_layeredgraph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_layeredgraph.cpp -I$(INCLUDE) -o $(OBJ)/gft_layeredgraph.o

$(OBJ)/gft_heap.o: $(SRC)/gft_heap.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap.o

$(OBJ)/gft_heap32.o: $(SRC)/gft_heap32.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap32.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap32.o

$(OBJ)/gft_heap64f.o: $(SRC)/gft_heap64f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap64f.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap64f.o

$(OBJ)/gft_heap_lex.o: $(SRC)/gft_heap_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap_lex.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap_lex.o

$(OBJ)/gft_heap32fi_lex.o: $(SRC)/gft_heap32fi_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap32fi_lex.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap32fi_lex.o

$(OBJ)/gft_heap32fif_lex.o: $(SRC)/gft_heap32fif_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap32fif_lex.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap32fif_lex.o

$(OBJ)/gft_heap32fifi_lex.o: $(SRC)/gft_heap32fifi_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap32fifi_lex.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap32fifi_lex.o

$(OBJ)/gft_heap32fiif_lex.o: $(SRC)/gft_heap32fiif_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap32fiif_lex.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap32fiif_lex.o

$(OBJ)/gft_heappair.o: $(SRC)/gft_heappair.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heappair.cpp -I$(INCLUDE) -o $(OBJ)/gft_heappair.o

$(OBJ)/gft_ift.o: $(SRC)/gft_ift.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_ift.cpp -I$(INCLUDE) -o $(OBJ)/gft_ift.o

$(OBJ)/gft_string.o: $(SRC)/gft_string.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_string.cpp -I$(INCLUDE) -o $(OBJ)/gft_string.o

$(OBJ)/gft_arraylist.o: $(SRC)/gft_arraylist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_arraylist.cpp -I$(INCLUDE) -o $(OBJ)/gft_arraylist.o

$(OBJ)/gft_filelist.o: $(SRC)/gft_filelist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_filelist.cpp -I$(INCLUDE) -o $(OBJ)/gft_filelist.o

$(OBJ)/gft_highlight.o: $(SRC)/gft_highlight.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_highlight.cpp -I$(INCLUDE) -o $(OBJ)/gft_highlight.o

$(OBJ)/gft_color.o: $(SRC)/gft_color.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_color.cpp -I$(INCLUDE) -o $(OBJ)/gft_color.o


$(OBJ)/gft_evaluation.o: $(SRC)/gft_evaluation.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_evaluation.cpp -I$(INCLUDE) -o $(OBJ)/gft_evaluation.o

$(OBJ)/gft_marker.o: $(SRC)/gft_marker.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_marker.cpp -I$(INCLUDE) -o $(OBJ)/gft_marker.o


$(OBJ)/gft_matrix.o: $(SRC)/gft_matrix.cpp $(INCLUDE)/gft_matrix.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_matrix.cpp -I$(INCLUDE) -o $(OBJ)/gft_matrix.o

$(OBJ)/gft_bmap.o: $(SRC)/gft_bmap.cpp $(INCLUDE)/gft_bmap.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_bmap.cpp -I$(INCLUDE) -o $(OBJ)/gft_bmap.o

$(OBJ)/gft_set.o: $(SRC)/gft_set.cpp $(INCLUDE)/gft_set.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_set.cpp -I$(INCLUDE) -o $(OBJ)/gft_set.o

$(OBJ)/gft_stack.o: $(SRC)/gft_stack.cpp $(INCLUDE)/gft_stack.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_stack.cpp -I$(INCLUDE) -o $(OBJ)/gft_stack.o

$(OBJ)/gft_analysis.o: $(SRC)/gft_analysis.cpp $(INCLUDE)/gft_analysis.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_analysis.cpp -I$(INCLUDE) -o $(OBJ)/gft_analysis.o

$(OBJ)/gft_pqueue32.o: $(SRC)/gft_pqueue32.cpp $(INCLUDE)/gft_pqueue32.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_pqueue32.cpp -I$(INCLUDE) -o $(OBJ)/gft_pqueue32.o

$(OBJ)/gft_pqueue16.o: $(SRC)/gft_pqueue16.cpp $(INCLUDE)/gft_pqueue16.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_pqueue16.cpp -I$(INCLUDE) -o $(OBJ)/gft_pqueue16.o

$(OBJ)/gft_hashtable.o: $(SRC)/gft_hashtable.cpp $(INCLUDE)/gft_hashtable.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_hashtable.cpp -I$(INCLUDE) -o $(OBJ)/gft_hashtable.o

$(OBJ)/gft_curve.o: $(SRC)/gft_curve.cpp $(INCLUDE)/gft_curve.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_curve.cpp -I$(INCLUDE) -o $(OBJ)/gft_curve.o

$(OBJ)/gft_radiometric.o: $(SRC)/gft_radiometric.cpp $(INCLUDE)/gft_radiometric.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_radiometric.cpp -I$(INCLUDE) -o $(OBJ)/gft_radiometric.o

$(OBJ)/gft_filtering.o: $(SRC)/gft_filtering.cpp $(INCLUDE)/gft_filtering.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_filtering.cpp -I$(INCLUDE) -o $(OBJ)/gft_filtering.o

$(OBJ)/gft_attributelist.o: $(SRC)/gft_attributelist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_attributelist.cpp -I$(INCLUDE) -o $(OBJ)/gft_attributelist.o

$(OBJ)/gft_configfile.o: $(SRC)/gft_configfile.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_configfile.cpp -I$(INCLUDE) -o $(OBJ)/gft_configfile.o

$(OBJ)/gft_geometry.o: $(SRC)/gft_geometry.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_geometry.cpp -I$(INCLUDE) -o $(OBJ)/gft_geometry.o


$(OBJ)/gft_draw.o: $(SRC)/gft_draw.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_draw.cpp -I$(INCLUDE) -o $(OBJ)/gft_draw.o


clean:
	rm -f $(LIB)/lib*.a; rm -f $(OBJ)/*.o $(BIN)/*


