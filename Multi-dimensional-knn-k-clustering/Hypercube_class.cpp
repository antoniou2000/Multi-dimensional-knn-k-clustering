#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "Hypercube_class.h"
#include "functions.h"
#include "data_handling.h"

using namespace std;

Hypercube::Hypercube(vector<pointType>& p, int k, int d, int w, int m, int probes, metricFunction f): k(k), d(d), w(w), probes(probes), M(m), points(p), metric(f){

    hashTable = new HashTable_Hypercube(k, d, w);

    for (int i = 0; i < p.size(); i++){
        hashTable->insertPoint(&p[i]);
    }

}

Hypercube::~Hypercube(){

    delete hashTable;

}
        
void Hypercube::nearestNeighbours(pointType* query, int N, resultNN& rNN){

    int counter = 0, query_id;
    double dist;
    bool flag = false;

    vector<pointType*> bucket;

    vector<double> dTrue(N); //True best distance
    vector<pair<double, string>> bestDist(N); //for each q {best distance, best candidate}

    //hypercube vertices that will be sorted using their hamming distance-> First place is Hamming Distance, Second is index
    vector<pair<int, int>> vertices(exp2(k));

    //Initialize values for each loop
    for (int i = 0; i < bestDist.size(); i++){
        bestDist[i] = {numeric_limits<double>::max(), "EMPTY"};
        dTrue[i] = numeric_limits<double>::max();
    }

    auto startTrue=chrono::high_resolution_clock::now();

    //Calculate true best distances (brute force) and store the distances (Not the points)
    double distanceTrue;
    for(int point = 0; point < points.size(); point++){
        distanceTrue = metric(&points[point], query, d);
        dTrue.insert(upper_bound(dTrue.cbegin(),dTrue.cend(),distanceTrue),distanceTrue);
        dTrue.pop_back();
    }

    auto stopTrue = chrono::high_resolution_clock::now();
    auto durationTrue = chrono::duration_cast<chrono::microseconds>(stopTrue-startTrue);

    auto start = chrono::high_resolution_clock::now();

    hashTable->getBucket(query, bucket); //get bucket to initialize query's hash_id_hypercube
    query_id = query->hash_id_hypercube;

    //Find hamming distances between query and each bucket and sort them
    for(int i=0;i<vertices.size();i++){
        vertices[i].first=hammingDistance(i,query_id);
        vertices[i].second=i;
    }
    sort(vertices.begin(),vertices.end());

    //Check |probe| vertices (+ the one in which q "belongs")
    for (int i = 0; i < probes; i++){

        if (counter == M){
            break;
        }

        flag = false;

        hashTable->getBucket(vertices[i].second, bucket);

        for (int j = 0; j < bucket.size(); j++){

            if (bucket[j]->hash_id_hypercube == query_id){

                counter++;
                flag = true;

                dist = metric(bucket[j], query, d);

                for(int index = 0; index < N; index++){
                        
                    //If point already in the bucket, ignore it
                    if (bucket[j]->id == bestDist[index].second){
                        break;
                    }

                    if(dist<bestDist[index].first){
                        //insert sorted in bucket {distance from point_j, name of point_j}
                        bestDist.insert(upper_bound(bestDist.begin(), bestDist.end(), make_pair(dist, bucket[j]->id)), make_pair(dist, bucket[j]->id));
                        bestDist.pop_back();
                        break;
                    }
                }

            }

            if (counter == M){
                break;
            }

        }

        //If no such points where found (*), then calculate the distance of each point in the bucket
        if (flag == false){

            for(int j=0;j<(int)bucket.size();j++){

                counter++;

                dist=metric(bucket[j], query, d);
                
                for(int index = 0; index < N; index++){

                    //If point already in the bucket, ignore it
                    if (bucket[j]->id == bestDist[index].second){
                        break;
                    }

                    //insert sorted in bucket {distance from point_j, name of point_j}
                    if(dist<bestDist[index].first){
                        bestDist.insert(upper_bound(bestDist.begin(), bestDist.end(), make_pair(dist, bucket[j]->id)), make_pair(dist, bucket[j]->id));
                        bestDist.pop_back();
                        break;
                    }
                }

                if (counter == M){
                    break;
                }

            } 

        }

    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop-start);

    rNN.Dist = bestDist;
    rNN.trueDist = dTrue;
    rNN.tDist = duration;
    rNN.tTrue = durationTrue;

}

void Hypercube::rangeSearch(pointType* query, int R, resultR& rR){

    int counter = 0, query_id;
    double dist;
    bool flag = false;

    vector<pointType*> bucket;

    vector<pointType*> dR; //points within r distance from each q

    //hypercube vertices that will be sorted using their hamming distance-> First place is Hamming Distance, Second is index
    vector<pair<int, int>> vertices(exp2(k));

    hashTable->getBucket(query, bucket);
    query_id = query->hash_id_hypercube;

    for(int i=0;i<vertices.size();i++){
        vertices[i].first=hammingDistance(i,query_id);
        vertices[i].second=i;
    }
    sort(vertices.begin(),vertices.end());

    for (int i = 0; i < probes; i++){

        if (counter == M){
            break;
        }

        hashTable->getBucket(vertices[i].second, bucket);

        for(int j = 0; j < bucket.size(); j++){

            counter++;

            dist=metric(bucket[j], query, d);

            //If distance smaller than R, insert it in dR vector
            if(dist <= R){

                //If element already in vector, skip
                if (find(dR.begin(), dR.end(), bucket[j]) != dR.end()){
                    continue;
                }

                //dR.insert(upper_bound(dR.cbegin(),dR.cend(),bucket[j]->id),bucket[j]->id);
                dR.push_back(bucket[j]);
            }

            if (counter == M){
                break;
            }

        }

    }

    rR.RDist = dR;

}
