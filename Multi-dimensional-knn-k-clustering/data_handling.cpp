#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

#include "data_handling.h"

std::vector<pointType> get_data(std::string fileName){

    int d = getNumOfDimensionsInFile(fileName);
    int numOfPoints = getNumOfLines(fileName);

    std::vector<pointType> temp(numOfPoints);
    for (int i = 0; i < numOfPoints; i++){
        temp[i].coords = new double[d]; //each point has d coords -> allocate memory
    }

    std::ifstream data;
    data.open(fileName);    //open file
    while(!data){
        std::cout<<"Could not open file. Please try again\nEnter your input file name\n-->";
        std::cin>>fileName;
        data.open(fileName); //open file
    }
    
    int i = 0; //i: the point that is currently being read - the ith line of the file
    std::string coords;
    int row=0;
    while (data >> coords ){ //while there are numbers in the file

        if(row==0){ //1st element of each row holds the ID/name of the data
            temp[i].id = coords;
        }
        else{ //Store the coordinate to our array
            temp[i].coords[row-1] = std::stod(coords);
        }

        if(row==d){ //finished reading the line, go to the next one and start over
            i++; //next line/point
            row=0; //reset the row counter
        }
        else{
            row++;
        }
    }
    data.close(); //close the file    
    std::cout<<"Data succesfully stored\n";
    return temp;
}

//items per line minus 1 (aka id)
int getNumOfDimensionsInFile(std::string File){
    int numOfDims=0;
    std::ifstream data;
    std::string line;
    data.open(File);
    getline(data,line);
    std::stringstream stream(line);
    while(stream>>line){
        numOfDims++;
    }
    return numOfDims-1;
}

//number of lines of file
int getNumOfLines(std::string File){
    int counter=0;
    std::string line;
    std::ifstream file(File);
    while (getline(file, line)){
        counter++;
    }
    return counter;
}