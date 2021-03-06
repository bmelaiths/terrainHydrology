#include <stdio.h>

int main() {
    //gather inputs
    char i, j;
    fread(&i, 1, 1, stdin);
    fread(&j, 1, 1, stdin);

    //perform necessary computations
    int result = i << j;

    //export outputs
    fwrite(&result, 1, 1, stdout);

    return 0;
}