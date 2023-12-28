
INCLUDE=./include
SRC=./src
OBJ=./obj

all: HLOIFT_closeness runexp

#Compiladores
CC=gcc
CXX=g++

FLAGS= -O3 -Wall
#FLAGS= -Wall -g
LFLAGS = -lpthread -lz

#Bibliotecas
GFTLIB  = -L./lib/gft/lib -lgft
GFTFLAGS  = -I./lib/gft/include

#Rules
libgft:
	$(MAKE) -C ./lib/gft


$(OBJ)/image_io.o: $(SRC)/image_io.c
	$(CXX) $(FLAGS) $(GFTFLAGS) -c $(SRC)/image_io.c -I$(INCLUDE) -o $(OBJ)/image_io.o

$(OBJ)/hloift.o: $(SRC)/hloift.c
	$(CXX) $(FLAGS) $(GFTFLAGS) -c $(SRC)/hloift.c -I$(INCLUDE) -o $(OBJ)/hloift.o

HLOIFT_closeness: \
HLOIFT_closeness.c \
$(OBJ)/image_io.o \
$(OBJ)/hloift.o libgft
	$(CXX) $(FLAGS) $(GFTFLAGS) -I$(INCLUDE) HLOIFT_closeness.c $(OBJ)/*.o \
	$(GFTLIB) -o HLOIFT_closeness -lm $(LFLAGS)

runexp: \
runexp.c libgft
	$(CXX) $(FLAGS) $(GFTFLAGS) -I$(INCLUDE) runexp.c $(OBJ)/*.o \
	$(GFTLIB) -o runexp -lm $(LFLAGS)

clean:
	$(RM) *~ *.o $(OBJ)/*.o HLOIFT_closeness runexp -f
