#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>

#include "LSH_class.h"
#include "functions.h"
#include "data_handling.h"

using namespace std;

LSH::LSH(vector<pointType>& p, int k, int l, int d, int w, metricFunction m): k(k), L(l), d(d), w(w), points(p), metric(m){

    //Allocate dynamically memory for the points to store the id(p) value
    // (for each different hash table aka for each different g() function, there is a different 
    // id() function, store them so as to not compute them constantly)
    // *This id(p) is stored while calling the HashTable_LSH::insertPoint function
    for (int i = 0; i < p.size(); i++){ 
        p[i].hash_id = new long[L];
    }

    //Hash table's size: n/8 or n/16 -> n/32
    int tableSize = p.size()/32;

    hashTables = new HashTable_LSH*[L];

    //Initialize L hash tables, and insert each point in every one of them
    for (int i = 0; i < L; i++){
        
        hashTables[i] = new HashTable_LSH(tableSize, k, d, w);

        for (int j = 0; j < p.size(); j++){
            hashTables[i]->insertPoint(&p[j]);
        }

    }

}

//Free all memory that was previously allocated
LSH::~LSH(){

    for (int i = 0; i < points.size(); i++){
        delete[] points[i].hash_id;
    }

    for (int i = 0; i < L; i++){
        delete hashTables[i];
    }
    delete[] hashTables;

}

//void LSH::nearestNeighbours(pointType query, resultNN& rNN){
void LSH::nearestNeighbours(pointType* query, int N, resultNN& rNN){

    double dist, distTrue;
    int query_id, numOfHashTables = L;
    bool flag = false;

    vector<pointType*> bucket;
    vector<double> dTrue(N); //True best distance
    vector<pair<double, string>> bestDist(N); //for each q {best distance, best candidate}
    
    //Initialiazation
    for (int i = 0; i < N; i++){
        bestDist[i] = {numeric_limits<double>::max(), "EMPTY"};
        dTrue[i] = numeric_limits<double>::max();
    }

    auto startTrue = chrono::high_resolution_clock::now();

    //Calculate true best distances (brute force) and store the distances (Not the points)
    for (int p = 0; p < points.size(); p++){
        distTrue = metric(&points[p], query, d);
        dTrue.insert(upper_bound(dTrue.cbegin(),dTrue.cend(),distTrue),distTrue);
        dTrue.pop_back();
    }

    auto stopTrue=chrono::high_resolution_clock::now();
    auto durationTrue=chrono::duration_cast<chrono::microseconds>(stopTrue-startTrue);

    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < numOfHashTables; i++){

        hashTables[i]->getBucket(query, bucket);
        query_id = query->hash_id[i];

        flag = false;

        for (int j = 0; j < bucket.size(); j++){

            //Initially calculate distances from points in the bucket with same id from q (*)
            if (bucket[j]->hash_id[i] == query_id){

                flag = true;
                dist = metric(bucket[j], query, d);

                for (int index = 0; index < N; index++){

                    //If point already in the bucket, ignore it
                    if (bucket[j]->id == bestDist[index].second){
                        break;
                    }

                    //insert sorted in bucket {distance from point_j, name of point_j}
                    if (dist < bestDist[index].first){
                        bestDist.insert(upper_bound(bestDist.begin(), bestDist.end(), make_pair(dist, bucket[j]->id)), make_pair(dist, bucket[j]->id));
                        bestDist.pop_back();
                        break;
                    }

                }

            }

        }

        //If no such points where found (*), then calculate the distance of each point in the bucket
        if (flag == false){

            for (int j = 0; j < bucket.size(); j++){

                dist = metric(bucket[j], query, d);

                for (int index = 0; index < N; index++){

                    //If point already in the bucket, ignore it
                    if (bucket[j]->id == bestDist[index].second){
                        break;
                    }

                    //insert sorted in bucket {distance from point_j, name of point_j}
                    if (dist < bestDist[index].first){
                        bestDist.insert(upper_bound(bestDist.begin(), bestDist.end(), make_pair(dist, bucket[j]->id)), make_pair(dist, bucket[j]->id));
                        bestDist.pop_back();
                        break;
                    }

                }

            }

        }

                

    }

    auto stop=chrono::high_resolution_clock::now();
    auto duration=chrono::duration_cast<chrono::microseconds>(stop-start);

    rNN.Dist = bestDist;
    rNN.trueDist = dTrue;
    rNN.tDist = duration;
    rNN.tTrue = durationTrue;

}

void LSH::rangeSearch(pointType* query, int R, resultR& rR){
    
    double dist;
    int query_id, numOfHashTables = L;

    vector<pointType*> bucket;
    vector<pointType*> dR; //points within r distance from each q

    for (int i = 0; i < numOfHashTables; i++){

        hashTables[i]->getBucket(query, bucket);
        query_id = query->hash_id[i];

        for (int j = 0; j < bucket.size(); j++){

            dist = metric(bucket[j], query, d);

            //If distance smaller than R, insert it in dR vector
            if (dist <= R){

                //If element already in vector, skip
                if (find(dR.begin(), dR.end(), bucket[j]) != dR.end()){
                    continue;
                }

                dR.push_back(bucket[j]);
            }

        }

    }

    rR.RDist = dR;
}