CC = g++
CFLAGS = -g -Wall -I.

OBJ_LSH = main_LSH.o HashTable_LSH.o data_handling.o functions.o LSH_class.o
OBJ_HC = main_Hypercube.o HashTable_Hypercube.o Hypercube_class.o functions.o data_handling.o
OBJ_CL = main_Clustering.o clustering.o LSH_class.o Hypercube_class.o HashTable_LSH.o HashTable_Hypercube.o functions.o data_handling.o

lsh: $(OBJ_LSH)
	$(CC) -o $@ $^ $(CFLAGS)

cube: $(OBJ_HC)
	$(CC) -o $@ $^ $(CFLAGS)

cluster: $(OBJ_CL)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o lsh cube cluster