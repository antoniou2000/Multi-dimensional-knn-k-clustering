#ifndef DATA_HANDLING_H

#define DATA_HANDLING_H

#include <vector>
#include "types.h"

std::vector<pointType> get_data(std::string);
int getNumOfDimensionsInFile(std::string);
int getNumOfLines(std::string);

#endif