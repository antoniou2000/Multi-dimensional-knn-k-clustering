#ifndef TYPES_H

#define TYPES_H

#include <string>
#include <vector>
#include <chrono>

#define M_VALUE 4294967291

typedef struct pointTag{
    std::string id;
    double *coords;
    long *hash_id;
    unsigned int hash_id_hypercube;
} pointType;

typedef struct hashParametersTag{
    int r;
    double *v;
    double t;
} hashParametersType;

typedef struct resultNNTag{
    std::vector<std::pair<double, std::string>> Dist;
    std::vector<double> trueDist;
    std::chrono::microseconds tDist;
    std::chrono::microseconds tTrue;
} resultNN;

typedef struct resultRTag{
    std::vector<pointType*> RDist;
} resultR;

typedef struct resultCLTag{
    std::vector<pointType> centroids;
    std::vector<std::vector<pointType*>> points_per_centroid;
    std::vector<double> silhouette;
    std::chrono::seconds time;
} resultCL;

#endif