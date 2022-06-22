#ifndef HASH_TABLE_H

#define HASH_TABLE_H

#include <vector>

#include "types.h"

//abstract Generic Hash Table with template type

template <typename T>
class GenericHashTable{

    private:

        std::vector<T> *table;

    protected:

        const int tableSize;
        virtual int hash_function(T) = 0;

    public:

        GenericHashTable(int s): tableSize(s){
            table = new std::vector<T>[tableSize];
        }
        ~GenericHashTable(){
            delete[] table;
        }

        void insertItem(T t){
            int pos = hash_function(t);
            table[pos].push_back(t);
        }

        void getBucket(T t, std::vector<T>& v){
            int pos = hash_function(t);
            v = table[pos];
        }

        void getBucket(int index, std::vector<T>& v){
            v = table[index];
        }
        
};

#endif