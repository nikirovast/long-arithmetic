//
// Created by Andrey Manucharyan on 20.06.21.
//
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

struct bignum {
    uint64_t size;
    uint64_t array[];
};

void sqrt2(uint64_t n, struct bignum* xn, struct bignum* xnp1);


int main(int argc, char** argv) {
    // implement parsing of arguments, size = 5 as example
    uint64_t size = 5;
    struct bignum* xn = malloc(size + 1); //size + 1 because of size parameter in the struct
    if (xn == NULL) {
        fprintf(stderr, "Not enough memory for the given number of iterations: %llu", size);
        return 1;
    }
    xn->size = 1;
    xn->array[0] = 1;

    struct bignum* xnp1 = malloc(size + 1); //size + 1 because of size parameter in the struct
    if (xnp1 == NULL) {
        fprintf(stderr, "Not enough memory for the given number of iterations: %llu", size);
        return 1;
    }
    xnp1->size = 1;
    xnp1->array[0] = 1;
    printf("memory for xn and xnp1 successfully allocated");
    return 0;
}

