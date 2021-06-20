//
// Created by Andrey Manucharyan on 20.06.21.
//
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "errno.h"

struct bignum {
    uint64_t size;
    uint64_t array[];
};

void sqrt2(uint64_t n, struct bignum* xn, struct bignum* xnp1);


int main(int argc, char** argv) {
    // implement parsing of arguments, size = 5 as example
    if (argc != 3) {
        fprintf(stderr, "usage: %s <number of iterations> <decimal (d) or hexadecimal (x) output>\n", argv[0]);
        return 1;
    }
    uint64_t size = strtoull(argv[1], NULL, 0);
    if (errno == ERANGE) {
        fprintf(stderr, "the given number %llu can not be represented, please pick a number < UINT64_MAX", size);
        return 1;
    }
    if (size == 0ULL) {
        fprintf(stderr, "invalid number of iterations: it should be > 0 or the given parameter was not a number");
        return 1;
    }
    char *output = argv[2];
    if ((*output != 'd' && *output != 'x') || strlen(output) > 1) {
        fprintf(stderr, "invalid numeral system for output");
        return 1;
    }
    struct bignum* xn = malloc(size + 1); //size + 1 because of size parameter in the struct
    if (xn == NULL) {
        fprintf(stderr, "not enough memory for the given number of iterations: %llu", size);
        return 1;
    }
    xn->size = 1;
    xn->array[0] = 1;

    struct bignum* xnp1 = malloc(size + 1); //size + 1 because of size parameter in the struct
    if (xnp1 == NULL) {
        fprintf(stderr, "not enough memory for the given number of iterations: %llu", size);
        return 1;
    }
    xnp1->size = 1;
    xnp1->array[0] = 2;
    printf("memory for xn and xnp1 successfully allocated, value in the first element of array with size %llu: xn is %llu and xnp1 is %llu.",size,  xn->array[0], xnp1->array[0]);
    return 0;
}

