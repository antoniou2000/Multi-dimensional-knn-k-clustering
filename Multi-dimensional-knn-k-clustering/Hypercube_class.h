#ifndef HYPERCUBE_CLASS_H

#define HYPERCUBE_CLASS_H

#include <vector>

#include "HashTable_Hypercube.h"
#include "types.h"

class Hypercube{

    //pointer to function for the metric function
    typedef double (*metricFunction)(pointType*, pointType*, int);

    private:
        std::vector<pointType>& points; //reference to points, so as to
            //free all the memory already allocated before finishing

        HashTable_Hypercube* hashTable;

        const int k; //k h() per hash table
        const int d; //d dimensions per point/query
        const int w; //w "window"
        const int probes;
        const int M;
        const metricFunction metric; //metric function

    public:
        Hypercube(std::vector<pointType>&, int, int, int, int, int, metricFunction);
        ~Hypercube();

        void nearestNeighbours(pointType*, int, resultNN&);
        void rangeSearch(pointType*, int, resultR&);

};

#endif