#ifndef CLUSTERING_H

#define CLUSTERING_H

#include <vector>
#include "types.h"
#include "LSH_class.h"
#include "data_handling.h"
#include "Hypercube_class.h"

class clustering{

    //pointer to function for the metric function
    typedef double (*metricFunction)(pointType*, pointType*, int);

    private:
        std::vector<pointType>& points; //reference to points, so as to
            //free all the memory already allocated before finishing

        std::vector<pointType> centroids;
        std::vector<std::vector<pointType*>> points_per_centroid;
        std::vector<double> silhouette;

        char method; // c -> classic (Lloyd's), l -> LSH, h -> Hypercube

        //parameters:
        const int K; // number of clusters/K-medians
        const int d; // number of dimensions per point

        int L; // number of hashtables (for LSH)

        int M; // max number of elements to be checked (for Hypercube)
        int probes; // max number of probes to be checked (for Hypercube)

        int k; // k number of h() function / f() (for LSH and Hypercube)

        metricFunction metric;

        //different algorithms:
        LSH *lsh;
        Hypercube *hc;

        std::vector<pointType> clusters;

        void initialize_common();
        void k_means();
        void assignPoint(pointType&, pointType);
        double minimum_distance(pointType, int&, int);

        void lloyds_algorithm();
        void rangeassignment_LSH();
        void rangeassignment_Hypercube();

    public:
        clustering(std::vector<pointType>&, int, int, metricFunction);
        clustering(std::vector<pointType>&, int, int, int, int, metricFunction);
        clustering(std::vector<pointType>&, int, int, int, int, int, metricFunction);
        ~clustering();

        void cluster(resultCL&);
        
};

#endif