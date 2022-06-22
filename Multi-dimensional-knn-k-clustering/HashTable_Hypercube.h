#ifndef HASH_TABLE_HYPERCUBE_H

#define HASH_TABLE_HYPERCUBE_H

#include "GenericHashTable.h"
#include "types.h"
#include <vector>
#include <map>

class HashTable_Hypercube : private GenericHashTable<pointType*>{

    private:
        int h(pointType*, int);

        int hash_function(pointType*);

        const int k;
        const int d;
        const int w;
        hashParametersType* h_p;

        //If an h(point) returns value h, then f(h) = 0/1 should be the
        // same for that h value
        //So every time a new h() value is encountered, search if there is already
        // the f(h) value stored in the map, else "flip a coin" and add it to the
        // map, so that further encounters of that value should have the same result
        //For each h() function a different map, so basically there are k maps
        std::map<int, bool>* f;

    public:

        HashTable_Hypercube(int, int, int);
        ~HashTable_Hypercube();

        void insertPoint(pointType*);
        void getBucket(pointType*, std::vector<pointType*>&);
        void getBucket(int , std::vector<pointType*>& );
};

#endif

