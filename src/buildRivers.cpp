#include <stdio.h>

/*
A list of data structures that will be used:
* A KDTree that will be used to query node locations
* A vector that will store actual node structs (this is where
  the node information will be, such as elevation, priority,
  etc). The tree will merely store an index in the vector
* Some data structure to represent the shore
*/

int main() {
    //gather inputs
    char i, j;
    fread(&i, 1, 1, stdin);
    fread(&j, 1, 1, stdin);

    //perform necessary computations
    int result = i << j;

    //export outputs
    fwrite(&result, 1, 1, stdout);

    //free resources

    return 0;
}