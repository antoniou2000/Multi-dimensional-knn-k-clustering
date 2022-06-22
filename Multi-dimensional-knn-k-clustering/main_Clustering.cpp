#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

#include "clustering.h"
#include "data_handling.h"
#include "functions.h"

using namespace std;

int main(int argc, char *argv[]){

    //paramaters and their default values
    bool complete = false;
    int number_of_clusters;
    int L = 3, k_LSH = 4, M = 10, k_HC = 3, probes = 2;
    char method;
    string inputFile, configFile, outputFile;

    for (int i = 0; i < argc; i++){

        string argv_i = argv[i];
        if (argv_i == "-i"){
            inputFile = argv[i+1];
        }
        else if (argv_i == "-c"){
            configFile = argv[i+1];
        }
        else if (argv_i == "-o"){
            outputFile = argv[i+1];
        }
        else if (argv_i == "-complete"){
            complete = true;
        }
        else if (argv_i == "-m"){

            string next_argv = argv[i+1];
            if (next_argv == "Classic"){
                method = 'c';
            }
            else if (next_argv == "LSH"){
                method = 'l';
            }
            else{ // "Hypercube"
                method = 'h';
            }

        }

    }

    //get parameters config file
    ifstream data;
    data.open(configFile);    //open file
    while(!data){
        cout<<"Could not open file. Please try again\nEnter your input file name\n-->";
        cin>>configFile;
        data.open(configFile); //open file
    }
    string coords;
    while(data>>coords){
        if(coords=="number_of_clusters:"){
            data>>coords;
            number_of_clusters=stoi(coords);
        }
        else if(coords=="number_of_vector_hash_tables:"){
            data>>coords;
            L=stoi(coords);
        }
        else if(coords=="number_of_vector_hash_functions:"){
            data>>coords;
            k_LSH=stoi(coords);
        }
        else if(coords=="max_number_M_hypercube:"){
            data>>coords;
            M=stoi(coords);
        }
        else if(coords=="number_of_hypercube_dimensions:"){
            data>>coords;
            k_HC=stoi(coords);
        }
        else if(coords=="number_of_probes:"){
            data>>coords;
            probes=stoi(coords);
        }
        
    }

    int d = getNumOfDimensionsInFile(inputFile);
    vector<pointType> points = get_data(inputFile);

    clustering* cl;

    if (method == 'c')
        cl = new clustering(points, number_of_clusters, d, L2);
    else if (method == 'l')
        cl = new clustering(points, number_of_clusters, L, k_LSH, d, L2);
    else
        cl = new clustering(points, number_of_clusters, M, probes, k_HC, d, L2);

    resultCL rCL;

    cl->cluster(rCL);

    ofstream outfile(outputFile);

    outfile << "Algorithm: ";
    if (method == 'c')
        outfile << "Lloyds" << endl;
    else if (method == 'l')
        outfile << "Range Search LSH" << endl;
    else
        outfile << "Range Search Hypercube" << endl;
    
    for (int i = 0; i < number_of_clusters; i++){

        outfile << rCL.centroids[i].id;
        outfile << "{size: " << rCL.points_per_centroid[i].size() << ", ";

        for (int c = 0; c < d; c++){
            outfile << rCL.centroids[i].coords[c] << " ";
        }
        outfile << "}" << endl;

    }
    outfile << "clustering_time: " << rCL.time.count() << " seconds" << endl;
    outfile << "Silhouette: [";
    for (int i = 0; i < number_of_clusters; i++){
        outfile << rCL.silhouette[i] << ", ";
    }
    outfile << rCL.silhouette[number_of_clusters] << "]" << endl;

    outfile << endl;
    if (complete == true){

        for (int i = 0; i < number_of_clusters; i++){

            outfile << rCL.centroids[i].id << " {";
            for (int c = 0; c < d; c++){
                outfile << rCL.centroids[i].coords[c] << " ";
            }

            for (int j = 0; j < rCL.points_per_centroid[i].size(); j++){
                outfile << ", " << rCL.points_per_centroid[i][j]->id;
            }
            outfile << "}" << endl;

        }

        outfile << endl;

    }

    delete cl;
    

   

}