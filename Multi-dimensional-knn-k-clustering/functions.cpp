#include <iostream>
#include <vector>
#include <cmath>

#include "functions.h"

using namespace std;

//extend mod for negative numbers
long mod(long D, long d){

    if (D < 0){
        return (-D)%d;
    }
    else{
        return D%d;
    }

}

//L2 metric
double L2(pointType* p1, pointType* p2, int dims){
    double distance=0;
    for(int i=0;i<dims;i++){
        distance +=(p1->coords[i]-p2->coords[i])*(p1->coords[i]-p2->coords[i]);
    }
    return sqrt(distance);
}

//hamming distance between two numbers
int hammingDistance(int n1, int n2){
    int x = n1 ^ n2;
    int setBits = 0;
 
    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }
 
    return setBits;
}