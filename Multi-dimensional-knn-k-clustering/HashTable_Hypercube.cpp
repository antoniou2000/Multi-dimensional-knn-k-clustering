#include <iostream>
#include <vector>
#include <random>

#include "HashTable_Hypercube.h"
#include "functions.h"

using namespace std;

HashTable_Hypercube::HashTable_Hypercube(int k, int d, int w): GenericHashTable(int(exp2(k))), k(k), d(d), w(w){

    //Initialization of the h() functions' parameters
    //For each of the h1, h2, ..., hk functions, we have to store
    // its parameters, aka the value t and the vector v, as well as
    // the value r_i, with which we multiply each h so as to get
    // the value of the g() function
    h_p = new hashParametersType[k];

    random_device r;

    default_random_engine gen(r());
    normal_distribution<double> n(0.0, 1.0);
    uniform_int_distribution<int> uni_i(0, numeric_limits<int>::max());
    uniform_real_distribution<double> uni_d(0.0, w);

    for (int i = 0; i < k; i++){

        h_p[i].v = new double[d];
        for (int j = 0; j < d; j++){
            h_p[i].v[j] = n(gen);
        }

        h_p[i].r = uni_i(gen);
        h_p[i].t = uni_d(gen);

    }

    f = new map<int, bool>[k];

}

HashTable_Hypercube::~HashTable_Hypercube(){

    for (int i = 0; i < k; i++){
        delete[] h_p[i].v;
    }
    delete[] h_p;
    delete[] f;

}

void HashTable_Hypercube::insertPoint(pointType* p){

    p->hash_id_hypercube = hash_function(p);
    GenericHashTable::insertItem(p);

}

void HashTable_Hypercube::getBucket(pointType* q, vector<pointType*>& v){

    q->hash_id_hypercube = hash_function(q);
    GenericHashTable::getBucket(q, v);

}

void HashTable_Hypercube::getBucket(int index, vector<pointType*>& v){
    GenericHashTable::getBucket(index, v);
}

int HashTable_Hypercube::hash_function(pointType* p){

    int current;
    bool f_value;
    map<int, bool>::iterator it;
    unsigned long concat = 0;
    random_device r;

    default_random_engine gen(r());
    uniform_int_distribution<int> uni_f(0, 1);
    for (int i = 0; i < k; i++){
        concat = concat << 1;

        current = h(p, i);

        it = f[i].find(current);
        //f value for current has already been calculated previously
        // by another point so use that one
        if (it != f[i].end()){ 
            f_value = it->second;
        }
        //we haven't encountered current previously so randomize f and
        // store it in the map
        else{
            f_value = uni_f(gen);
            f[i][current] = f_value;
        }

        concat = concat | f_value;
        
    }
    
    return concat;

}

//h_i function for point p
int HashTable_Hypercube::h(pointType* p, int h_index){

    double p_v = 0.0;

    for (int i = 0; i < d; i++){
        p_v += (p->coords[i] * h_p[h_index].v[i])/w;
    }

    return (int)(p_v + (h_p[h_index].t / w));
}

