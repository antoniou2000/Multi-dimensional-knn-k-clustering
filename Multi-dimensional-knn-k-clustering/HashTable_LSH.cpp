#include <iostream>
#include <vector>
#include <random>
#include <climits>

#include "HashTable_LSH.h"
#include "functions.h"

using namespace std;

int HashTable_LSH::count{0}; //Initialize the counter

HashTable_LSH::HashTable_LSH(int s, int k, int d, int w): GenericHashTable(s), k(k), d(d), w(w){

    //Initialization of the h() functions' parameters
    //For each of the h1, h2, ..., hk functions, we have to store
    // its parameters, aka the value t and the vector v, as well as
    // the value r_i, with which we multiply each h so as to get
    // the value of the g() function
    h_p = new hashParametersType[k];

    random_device r;

    default_random_engine gen(r());
    normal_distribution<double> n(0.0, 1.0);
    uniform_int_distribution<int> uni_i(0, INT_MAX);
    uniform_real_distribution<double> uni_d(0.0, w);

    for (int i = 0; i < k; i++){

        h_p[i].v = new double[d];
        for (int j = 0; j < d; j++){
            h_p[i].v[j] = n(gen);
        }

        h_p[i].r = uni_i(gen);
        h_p[i].t = uni_d(gen);

    }

    tag = count++; //get the hash table's "id"

}

HashTable_LSH::~HashTable_LSH(){

    for (int i = 0; i < k; i++){
        delete[] h_p[i].v;
    }
    delete[] h_p;

}

void HashTable_LSH::insertPoint(pointType* p){

    p->hash_id[tag] = id(p); //store the id so as to not recalculate
        //every time

    GenericHashTable::insertItem(p);

}

void HashTable_LSH::getBucket(pointType* q, vector<pointType*>& v){

    q->hash_id[tag] = id(q); //for query q, calculate its id and
        // store it to avoid recalculating in future uses
    GenericHashTable::getBucket(q, v);

}

//g() function
int HashTable_LSH::hash_function(pointType* p){

    return mod(p->hash_id[tag], tableSize);

}

//ID(p) function
long HashTable_LSH::id(pointType* p){

    unsigned long id_value = 0;
    unsigned int h_i, a, b;

    for (int i = 0; i < k; i++){
        h_i = h(p, i);
        
        a = mod(h_i, M_VALUE);
        b = mod(h_p[i].r, M_VALUE);

        id_value += mod((a + b),M_VALUE);
    }

    return id_value;
}

//h_i function for point p
int HashTable_LSH::h(pointType* p, int h_index){

    double p_v = 0.0;

    for (int i = 0; i < d; i++){
        p_v += (p->coords[i] * h_p[h_index].v[i])/w;
    }

    return (int)(p_v + (h_p[h_index].t / w));
}