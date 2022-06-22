#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

#include "LSH_class.h"
#include "data_handling.h"
#include "functions.h"

using namespace std;

int main(int argc, char *argv[]){

    bool inputHasValue = false, queryHasValue = false, outputHasValue = false;
    string inputFile, queryFile, outputFile;

    int k = 4, L = 5, N = 1;
    double R = 500;
    char c;

    //Read from command line
    for (int i = 0; i < argc; i++){

        string argv_i = argv[i];
        if (argv_i == "-i"){
            inputHasValue = true;
            inputFile = argv[i+1];
        }
        else if (argv_i == "-q"){
            queryHasValue = true;
            queryFile = argv[i+1];
        }
        else if (argv_i == "-k"){
            k = atoi(argv[i+1]);
        }
        else if (argv_i == "-L"){
            L = atoi(argv[i+1]);
        }
        else if (argv_i == "-o"){
            outputHasValue = true;
            outputFile = argv[i+1];
        }
        else if (argv_i == "-N"){
            N = atoi(argv[i+1]);
        }
        else if (argv_i == "-R"){
            R = atoi(argv[i+1]);
        }
    }

    if (!inputHasValue){
        cout << "Enter path for dataset: ";
        cin >> inputFile;
        inputHasValue = true;
    }

    int d = getNumOfDimensionsInFile(inputFile);
    vector<pointType> points = get_data(inputFile);
    LSH lsh(points, k, L, d, 6, L2);

    vector<pointType> queries;

    while(1){

        if (!queryHasValue){
            cout << "Enter path for query file: ";
            cin >> queryFile;
            queryHasValue = true;
        }

        if (!outputHasValue){
            cout << "Enter path for output file: ";
            cin >> outputFile;
            outputHasValue = true;
        }

        ofstream outfile(outputFile);

        vector<pointType> queries;
        queries = get_data(queryFile);
        for (int q = 0; q < queries.size(); q++){ 
            queries[q].hash_id = new long[L];
        }

        for (int q = 0; q < queries.size(); q++){

            cout << "Searching for nearest neighbors for query " << queries[q].id << "...";

            resultNN rNN;
            resultR rR;

            lsh.nearestNeighbours(&queries[q], N, rNN);
            lsh.rangeSearch(&queries[q], R, rR);

            outfile << "Query: " << queries[q].id << endl;
            for(int i=0;i<N;i++){
                outfile << "Nearest Neighbor-" <<i<< ": ";
                outfile << rNN.Dist[i].second << endl;
                outfile << "distanceLSH: "<< rNN.Dist[i].first << endl;
                outfile << "distanceTrue: "<< rNN.trueDist[i] << endl;
            }
            outfile << "tLSH: " << rNN.tDist.count() << endl;
            outfile << "tTrue: " << rNN.tTrue.count() << endl;
            outfile << "R-near neighbours" << endl;
            for (int i = 0; i < rR.RDist.size(); i++){
                outfile << rR.RDist[i]->id << endl;
            }
            outfile << endl;

            cout << "Done" << endl;

        }

        outfile.close();

        for(int q = 0; q < queries.size(); q++){
            delete[] queries[q].hash_id;
            delete[] queries[q].coords;
        }

        cout << "Want to continue for another query file? (y/n): ";
        cin >> c;

        if (c == 'n')
            break;

        else{
            queryHasValue = false;
            outputHasValue = false;
        }

    }
    
    for (int i = 0; i < points.size(); i++){
        delete[] points[i].coords;
    }


}
