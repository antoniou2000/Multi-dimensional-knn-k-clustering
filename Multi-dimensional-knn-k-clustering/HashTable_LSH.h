#ifndef HASH_TABLE_LSH_H

#define HASH_TABLE_LSH_H

#include "GenericHashTable.h"
#include "types.h"
#include <vector>

class HashTable_LSH : private GenericHashTable<pointType*>{

    private:
        int tag; //The "id" of the hash table (0-L)
            // used to access, for each point p, its correct ID(p) function
            // (ID(p) = r1*h1 + .. rk*hk mod M)
        static int count; //for this reason, keep a static counter to
            // get the correct tag for this hash table

        //h_i function, using the values v_i, t_i (different for each h)
        // so as to compute the main hash function, g()
        int h(pointType*, int);

        int hash_function(pointType*);

        const int k; //k h() functions per g()
        const int d;
        const int w;
        hashParametersType* h_p; //and for each h(), keep its parameters

    public:

        HashTable_LSH(int, int, int, int);
        ~HashTable_LSH();

        void insertPoint(pointType*);
        void getBucket(pointType*, std::vector<pointType*>&);

        //id function: calculates the expression r1*h1 + ... + rk*hk mod M
        // it is used:
        //-once per point, when the point is first inserted in the hash table 
        // (then it's stored inside p.hash_id[tag] to avoid computing again)
        //-for the query q, when we're searching for its bucket
        long id(pointType*);

};

#endif