#include "errno.h"
#include "math.h"
#include <time.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TWO_POW_31 2147483648
#define N 32
#define ELEM_SIZE_MAX UINT32_MAX
#define LOG2_10 3.3
#define LOG10_2 0.3
#define POWER10_9 1000000000
typedef uint32_t elem_size_t;
typedef struct bignum bignum;
typedef struct matrix matrix;
static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};
struct bignum
{
    size_t size;
    elem_size_t *array;
};

/**
 * Calculate elapsed time
 */
double time_s(struct timespec s, struct timespec e) {
    return (double)(e.tv_sec - s.tv_sec) + (double) (e.tv_nsec - s.tv_nsec) * 1e-9;
}

extern void sum(uint64_t n, bignum *xn, bignum *xnp1);
void mulToPrint(elem_size_t *decArray, uint64_t *size);
bignum *mul(bignum *xn, bignum *xnp1);
void zeroJustify(bignum *n);
void arrayShift(bignum *n, uint64_t count);
int compareBignum(bignum *xn, bignum *xnp1);
bignum *div2bignums(bignum *xn, bignum *xnp1);
bignum *add(bignum *xn, bignum *xnp1);
bignum *sub(bignum *xn, bignum *xnp1);
bignum *mul2num(elem_size_t x, elem_size_t y);
matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2);
matrix *matrixBinaryExponentiation(uint64_t n, uint64_t highestBit);
uint64_t calculateHighestBit(uint64_t n);
uint64_t convertAccToN(uint64_t numDigits);
char *hexToPrint(bignum *a);
char *decToPrint(bignum *a);
void printUsage(char **argv);
void freeBigNum(bignum *toFree);
bignum *divideLongDivision(bignum *dividend, bignum *divisor);
bignum *copy(bignum *from, uint64_t countBits);
bignum *bitShiftRight(bignum *n, uint64_t count);
bignum *bitShiftLeft(bignum *n, uint64_t count);
bignum *addIntToBignum(bignum *n, uint64_t toAdd);
matrix *matrixSimpleExponentiation(uint64_t n);
void freeMatrix(matrix *toFree);
uint64_t countBits(elem_size_t *array, uint64_t *size);
char *toHex(bignum *a, uint64_t number);

bignum *new_bignum(size_t size) {
    bignum *res = malloc(sizeof(bignum));
    if (res == NULL) {
        fprintf(stderr, "Couldn't allocate memory for a struct\n");
        exit(1);
    }
    res->size = size;
    if (size == 0) {
        return res;
    }
    res->array = malloc(sizeof(elem_size_t) * size);
    if (res->array == NULL) {
        fprintf(stderr, "Couldn't allocate memory for an array\n");
        exit(1);
    }
    return res;
}

struct matrix {
    bignum *xn;
    bignum *xnp1;
    bignum *xnm1;
};

int main(int argc, char **argv) {
    int hexadecimal = 0;
    int c;
    int output = 0;
    if (argc == 1 || argc > 4) {
        printUsage(argv);
        return 1;
    }
    while (1) {
        int option_index = 0;
        c = getopt_long(argc, argv, "hxt", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
            case 'h': {
                printUsage(argv);
                return 1;
            }
            case 'x': {
                hexadecimal = 1;
                break;
            }
            case 't': {
                output = 1;
                break;
            }
            default: {
                fprintf(stderr, "Unknown option was used\n");
                exit(1);
            }
        }
    }
    // + flags because of a known bug of getops_long changing the order of arguments
    uint64_t n = strtoull(argv[1 + hexadecimal + output], NULL, 0);
    //we'll need to add 3 to n later if hexadecimal output
    if (errno == ERANGE || (hexadecimal && n >= UINT64_MAX - 2)) {
        fprintf(stderr, "the given number can not be represented, please pick a number <= %lu\n", UINT64_MAX);
        return 1;
    }
    if (n == 0ULL || *argv[1 + hexadecimal + output] == '-') {
        fprintf(stderr, "invalid number of digits: it should be > 0 or the given parameter was not a number\n");
        return 1;
    }
    if (*argv[1] == '0') {
        fprintf(stderr, "Given number begins with 0, enter a number 0 < number <= %lu\n", UINT64_MAX);
        return 1;
    }

    uint64_t op = convertAccToN(n);
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);
    matrix *res = matrixBinaryExponentiation(op, calculateHighestBit(op));
    //one might use the next line as comparison to the binary exponentiation
    //matrix *res = matrixSimpleExponentiation(op);
    if (output && hexadecimal) {
        char *xn = hexToPrint(res->xn);
        char *xnp1 = hexToPrint(res->xnp1);
        printf("x_%lu = %s\nx_%lu = %s\n", n, xn, n + 1, xnp1);
        free(xn);
        free(xnp1);
    }
    else if (output) {
        char *xn = decToPrint(res->xn);
        char *xnp1 = decToPrint(res->xnp1);
        printf("x_%lu = %s\nx_%lu = %s\n", n, xn, n + 1, xnp1);
        free(xn);
        free(xnp1);
    }
    uint64_t input = n;
    // we need it in order to convert decimal to hexadecimal which is pretty complex
    if (hexadecimal) {
        n += 3;
    }
    struct timespec end;
    bignum *oper = new_bignum(1);
    *oper->array = 10;
    for (uint64_t i = 0; i < n; i++) {
        bignum *savePtr = mul(res->xn, oper);
        freeBigNum(res->xn);
        res->xn = savePtr;
    }
    bignum *div = div2bignums(res->xn, res->xnp1);
    clock_gettime(CLOCK_MONOTONIC, &end);
    char *final;
    if (hexadecimal) {
        final = toHex(div, n - 3);
    } else {
        final = decToPrint(div);
    }
    printf("Took %f seconds to calculate x_%lu, x_%lu and divide them\n", time_s(start, end), input, input + 1);
    printf("Result after division: 1,%s\n", final);
    free(final);
    freeBigNum(div);
    freeBigNum(oper);
    freeBigNum(res->xnm1);
    freeBigNum(res->xnp1);
    free(res);
    //freeMatrix(res);
    return 0;
}


void printUsage(char **argv)
{
    printf("Usage: %s <number of digits> \n", argv[0]);
    printf("-x: sets output numeral system to hexadecimal, default: decimal\n");
    printf("-t: two numbers xn and xnp1 will appear on the screen the division of which results in sqrt(2) - 1\n");
    printf("-h|--help: for usage help\n");
}

/**
*
* Multiply two matrices.
*
*/

matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2)
{
    matrix *res = malloc(3 * sizeof(struct bignum));
    if (res == 0) {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixMul\n");
        exit(1);
    }
    bignum *tmp1 = mul(matrix1->xnm1, matrix2->xnm1);
    bignum *tmp2 = mul(matrix1->xn, matrix2->xn);
    res->xnm1 = add(tmp1, tmp2);
    freeBigNum(tmp1);
    freeBigNum(tmp2);
    tmp1 = mul(matrix1->xnm1, matrix2->xn);
    tmp2 = mul(matrix1->xn, matrix2->xnp1);
    res->xn = add(tmp1, tmp2);
    freeBigNum(tmp1);
    freeBigNum(tmp2);
    tmp1 = mul(matrix1->xn, matrix2->xn);
    tmp2 = mul(matrix1->xnp1, matrix2->xnp1);
    res->xnp1 = add(tmp1, tmp2);
    freeBigNum(tmp1);
    freeBigNum(tmp2);
    return res;
}

matrix *matrixSimpleExponentiation(uint64_t n) {
    struct matrix *matrix = malloc(3 * sizeof(bignum));
    if (matrix == NULL) {
        fprintf(stderr, "Couldn't allocate memory for a matrix\n");
        exit(1);
    }
    matrix->xnp1 = new_bignum(1);
    matrix->xnp1->array[0] = 2;
    matrix->xn = new_bignum(1);
    matrix->xn->array[0] = 1;
    matrix->xnm1 = new_bignum(0);
    struct matrix *matrixInitial = malloc(3 * sizeof(bignum));
    if (matrixInitial == NULL) {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixBinaryExp\n");
        exit(1);
    }
    matrixInitial->xnm1 = new_bignum(0);
    matrixInitial->xn = new_bignum(1);
    *(matrixInitial->xn->array) = 1;
    matrixInitial->xnp1 = new_bignum(1);
    *(matrixInitial->xnp1->array) = 2;
    for (uint64_t i = 0; i < n - 1; i++) {
        struct matrix *tmp = matrixMultiplication(matrix, matrixInitial);
        freeMatrix(matrix);
        matrix = tmp;
    }
    freeBigNum(matrixInitial->xn);
    freeBigNum(matrixInitial->xnm1);
    freeBigNum(matrixInitial->xnp1);
    free(matrixInitial);
    return matrix;
}

/**
 *
 * Binary exponentiation of the matrix to the exponent n as described in pdf
 *
 */

matrix *matrixBinaryExponentiation(uint64_t n, uint64_t highestBit) {
    struct matrix *matrix = malloc(3 * sizeof(bignum));
    if (matrix == NULL) {
        fprintf(stderr, "Couldn't allocate memory for a matrix\n");
        exit(1);
    }
    matrix->xnp1 = new_bignum(1);
    matrix->xnp1->array[0] = 2;
    matrix->xn = new_bignum(1);
    matrix->xn->array[0] = 1;
    matrix->xnm1 = new_bignum(1);
    matrix->xnm1->array[0] = 0;
    struct matrix *matrixInitial = malloc(3 * sizeof(bignum));
    if (matrixInitial == NULL) {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixBinaryExp\n");
        exit(1);
    }
    matrixInitial->xnm1 = new_bignum(1);
    *(matrixInitial->xnm1->array) = 0;
    matrixInitial->xn = new_bignum(1);
    *(matrixInitial->xn->array) = 1;
    matrixInitial->xnp1 = new_bignum(1);
    *(matrixInitial->xnp1->array) = 2;
    uint64_t i = highestBit - 1;
    while (1) {
        struct matrix *tmp = matrixMultiplication(matrix, matrix);
        freeMatrix(matrix);
        matrix = tmp;
        if (n >> i & 1) {
            tmp = matrixMultiplication(matrix, matrixInitial);
            freeMatrix(matrix);
            matrix = tmp;
        }
        if (i == 0) {
            break;
        }
        i--;
    }
    freeBigNum(matrixInitial->xn);
    freeBigNum(matrixInitial->xnm1);
    freeBigNum(matrixInitial->xnp1);
    free(matrixInitial);
    return matrix;
}

/**
 *
 * Calculate highest bit of the number
 *
 */

uint64_t calculateHighestBit(uint64_t n) {
    int order = -1;
    for (int i = 0; i < N; i++) {
        if ((n >> i) & 1) {
            order = i;
        }
    }
    return order;
}

/**
 *
 * Sum two big numbers.
 *
 */

bignum *add(bignum *xn, bignum *xnp1) {
    int flag = 0;
    if (xn->size > xnp1->size) {
        bignum *tmp = xn;
        xn = xnp1;
        xnp1 = tmp;
        flag = 1;
    }

    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;
    bignum *res = new_bignum(size2);
    uint64_t i;
    elem_size_t carry = 0;
    for (i = 0; i < size1; i++) {
        elem_size_t a = *(xnp1->array + i);
        elem_size_t b = *(xn->array + i);
        *(res->array + i) = a + b + carry;
        carry = (a > ELEM_SIZE_MAX - b);
    }
    while (i < size2) {
        elem_size_t a = *(xnp1->array + i);
        *(res->array + i) = a + carry;
        carry = (a > ELEM_SIZE_MAX - carry);
        i++;
    }
    if (carry) {
        elem_size_t *tmp = realloc(res->array, sizeof(elem_size_t) * (size2 + 1));
        if (tmp == NULL) {
            fprintf(stderr, "Couldn't reallocate memory for a carry in sum\n");
            exit(1);
        } else {
            res->array = tmp;
            res->size += 1;
            *(res->array + i) = carry;
        }
    }
    if (flag != 0) {
        bignum *tmp = xn;
        xn = xnp1;
        xnp1 = tmp;
    }
    return res;
}

/**
 *
 * Subtract xnp1 from xn.
 *
 */

bignum *sub(bignum *xn, bignum *xnp1) {
    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;
    bignum *res = new_bignum(sizeof(elem_size_t) * size1);
    res->size = size1;
    uint64_t i;
    elem_size_t carry = 0;
    for (i = 0; i < size2; i++) {
        elem_size_t a = *(xn->array + i);
        elem_size_t b = *(xnp1->array + i);
        if ((a < 1 && (b > 0 || carry > 0)) || a - carry < b) {
            *(res->array + i) = ELEM_SIZE_MAX - b - carry + a + 1;
            carry = 1;
        } else {
            *(res->array + i) = a - b - carry;
            carry = 0;
        }
    }
    while (i < size1) {
        elem_size_t a = *(xn->array + i);
        if (a >= carry) {
            *(res->array + i) = a - carry;
            carry = 0;
        } else {
            *(res->array + i) = ELEM_SIZE_MAX - carry + a + 1;
            carry = 1;
        }
        i++;
    }
    zeroJustify(res);
    return res;
}

/**
 *
 * Multiply two big numbers xnp1 by xn.
 *
 */

bignum *mul(bignum *xn, bignum *xnp1) {
    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;

    if (xn->size == 0 || xnp1->size == 0 || (xn->size == 1 && *xn->array == 0) || (xnp1->size == 1 && *xnp1->array == 0))
    {
        bignum *res = malloc(sizeof(*xn));
        if (res == 0) {
            fprintf(stderr, "Couldn't allocate memory for a zero res in mul\n");
            exit(1);
        }
        res->size = 0;
        return res;
    } else if (size1 > 1 || size2 > 1) {
        uint64_t aNewSize;
        uint64_t bNewSize;
        uint64_t cNewSize;
        uint64_t dNewSize;
        if (size1 < size2) {
            aNewSize = size1 / 2;
            bNewSize = size1 - aNewSize;
            dNewSize = bNewSize;
            cNewSize = size2 - dNewSize;
        } else {
            cNewSize = size2 / 2;
            dNewSize = size2 - cNewSize;
            bNewSize = dNewSize;
            aNewSize = size1 - bNewSize;
        }
        bignum *b = malloc(sizeof(bignum));
        if (b == NULL) {
            fprintf(stderr, "Couldn't allocate memory for b in mul\n");
            exit(1);
        }
        b->size = bNewSize;
        // we need to write in b the second part of xn
        b->array = xn->array;

        bignum *a = malloc(sizeof(bignum));
        if (a == NULL) {
            fprintf(stderr, "Couldn't allocate memory for a in mul\n");
            exit(1);
        }
        a->size = aNewSize;
        a->array = xn->array + bNewSize;
        bignum *c = malloc(sizeof(bignum));
        if (c == NULL) {
            fprintf(stderr, "Couldn't allocate memory for c in mul\n");
            exit(1);
        }
        c->size = cNewSize;
        c->array = xnp1->array + dNewSize;

        bignum *d = malloc(sizeof(bignum));
        if (d == NULL) {
            fprintf(stderr, "Couldn't allocate memory for d in mul\n");
            exit(1);
        }
        d->size = dNewSize;
        d->array = xnp1->array;
        bignum *ac = mul(a, c);
        bignum *bd = mul(b, d);
        bignum *a_add_b = add(a, b);
        bignum *c_add_d = add(c, d);
        free(a);
        free(b);
        free(c);
        free(d);
        bignum *abcd = mul(a_add_b, c_add_d);
        freeBigNum(a_add_b);
        freeBigNum(c_add_d);
        bignum *abcd_sub_ac = sub(abcd, ac);
        bignum *ad_add_bc = sub(abcd_sub_ac, bd);
        if (ac->size != 0) {
            arrayShift(ac, bNewSize + dNewSize);
        }
        if (ad_add_bc->size != 0) {
            arrayShift(ad_add_bc, bNewSize);
        }
        bignum *savePtr = add(ac, ad_add_bc);
        bignum *res = add(savePtr, bd);
        freeBigNum(ac);
        freeBigNum(bd);
        freeBigNum(abcd);
        freeBigNum(abcd_sub_ac);
        freeBigNum(ad_add_bc);
        freeBigNum(savePtr);
        zeroJustify(res);

        return res;
    } else {
        return mul2num(*(xn->array), *(xnp1->array));
    }
}

/**
 *
 * Multiply two elem_size_t numbers x by y.
 *
 */

bignum *mul2num(elem_size_t x, elem_size_t y) {
    if (x == 0 || y == 0) {
        bignum *res = malloc(sizeof(*res));
        if (res == 0) {
            fprintf(stderr, "Couldn't allocate memory for a zero res in mul\n");
            exit(1);
        }
        res->size = 0;
        return res;
    } else {
        size_t shift_size = N / 2;
        elem_size_t a = x >> shift_size;
        elem_size_t b = (x << shift_size) >> shift_size;
        elem_size_t c = y >> shift_size;
        elem_size_t d = (y << shift_size) >> shift_size;
        bignum ac;
        ac.size = 2;
        elem_size_t ac_array[] = {0, a * c};

        ac.array = ac_array;

        bignum ad;
        ad.size = 2;
        elem_size_t ad_array[] = {((a * d) << shift_size), (a * d) >> shift_size};
        ad.array = ad_array;

        bignum bc;
        bc.size = 2;
        elem_size_t bc_array[] = {(b * c) << shift_size, (b * c) >> shift_size};
        bc.array = bc_array;

        bignum bd;
        bd.size = 1;
        elem_size_t bd_array[] = {b * d};
        bd.array = bd_array;
        bignum *ac_add_ad = add(&ac, &ad);
        bignum *ac_add_ad_add_bc = add(ac_add_ad, &bc);
        freeBigNum(ac_add_ad);
        bignum *res = add(ac_add_ad_add_bc, &bd);
        freeBigNum(ac_add_ad_add_bc);
        zeroJustify(res);
        return res;
    }
}

/**
 *
 * Divide big number xn by big number xnp1.
 *
 */

bignum *div2bignums(bignum *xn, bignum *xnp1) {
    bignum *res = new_bignum(1);
    res->array[0] = 0;
    uint64_t xnBits = xn->size >= 2 ? (xn->size - 1) * N + calculateHighestBit(xn->array[xn->size - 1]) + 1
                                    : calculateHighestBit(xn->array[0]) + 1;
    uint64_t xnp1Bits = xnp1->size >= 2 ? (xnp1->size - 1) * N + calculateHighestBit(xnp1->array[xnp1->size - 1] + 1)
                                        : calculateHighestBit(xnp1->array[0]) + 1;
    uint64_t diff = xnBits - xnp1Bits;
    while (1) {
        bignum *toAdd = new_bignum(1);
        toAdd->array[0] = 1;
        bitShiftLeft(toAdd, diff);
        bignum *temp2 = mul(xnp1, toAdd);
        int compare = compareBignum(xn, temp2);
        if (compare < 1) {
            bignum *savePtr = add(res, toAdd);
            freeBigNum(res);
            res = savePtr;
            savePtr = sub(xn, temp2);
            freeBigNum(xn);
            xn = savePtr;
            freeBigNum(temp2);
            freeBigNum(toAdd);
            if (compare == 0) {
                break;
            }
        } else {
            freeBigNum(temp2);
            freeBigNum(toAdd);
        }

        if (diff == 0) {
            break;
        }
        diff--;
    }
    freeBigNum(xn);
    return res;
}

int compareBignum(bignum *xn, bignum *xnp1) {
    if (xn->size > xnp1->size) {
        return -1;
    } else if (xnp1->size > xn->size) {
        return 1;
    }
    uint64_t i = xnp1->size - 1;
    while (1) {
        if (xn->array[i] > xnp1->array[i]) {
            return -1;
        }
        if (xn->array[i] < xnp1->array[i]) {
            return 1;
        }
        if (i == 0) {
            break;
        }
        i--;
    }
    return 0;
}

void arrayShift(bignum *n, uint64_t count) {
    uint64_t i;
    elem_size_t *tmp = realloc(n->array, (n->size + count) * sizeof(elem_size_t));
    if (tmp == NULL) {
        fprintf(stderr, "Couldn't reallocate memory for a n in arrayShift\n");
        exit(1);
    } else {
        n->array = tmp;
    }
    if ((n->size == 0) && (n->array[0] == 0))
        return;
    i = n->size - 1;
    while (1) {
        n->array[i + count] = n->array[i];
        if (i == 0) {
            break;
        }
        i--;
    }

    for (i = 0; i < count; i++) {
        *(n->array + i) = 0;
    }
    n->size = n->size + count;
    zeroJustify(n);
}

void zeroJustify(bignum *n) {
    while ((n->size == 1 && n->array[0] == 0) || ((n->size > 1) && (n->array[n->size - 1] == 0))) {
        n->size--;
    }
}

/**
 *
 * Counts num of iteratons needed to achieve certain accuracy
 *
 */

uint64_t convertAccToN(uint64_t numDigits) {
    return 2 + (LOG2_10 / 2) * numDigits;
}

/**
 *
 * Converts bignum to hexadecimal format string
 *
 */

char *hexToPrint(bignum *a) {
    uint64_t count = a->size;
    if (count == 0) {
        char *zero = malloc(2);
        if (zero == NULL) {
            fprintf(stderr, "Couldn't allocate memory for a zerostring\n");
            exit(1);
        }
        *zero = '0';
        *(zero + 1) = 0;
        return zero;
    }

    int numCharinOneElem = N / 4;
    count--;
    elem_size_t *array = a->array;
    // could be confusing: first in the decomal format, last elem of array
    elem_size_t first = *(array + count);
    // to check how many bytes we need for the last elem of array
    char *tmp = malloc(9);
    if (tmp == 0) {
        fprintf(stderr, "Couldn't allocate memory for a tmp\n");
        exit(1);
    }
    // always add plus 1 byte because snprintf always adds and of the string, but we overwrite it then
    uint16_t addLen = snprintf(tmp, 9, "%x", first);
    free(tmp);
    // One extra for a string terminator
    uint64_t len = count * numCharinOneElem + addLen + 1;

    // string a pointer to the beginning of the string,
    // str the pointer to the location where we should write right now
    char *string, *str;
    str = string = malloc(len);
    if (string == NULL) {
        fprintf(stderr, "Couldn't allocate memory for a string in hexPrint\n");
        exit(1);
    }
    snprintf(str, addLen + 1, "%x", first);
    str += addLen;
    if (count == 0) {
        *str = '\0';
        return string;
    }
    count--;
    while (1) {
        elem_size_t c = *(array + count);
        snprintf(str, numCharinOneElem + 1, "%08x", c);
        str += numCharinOneElem;
        if (count == 0) {
            break;
        } else {
            count--;
        }
    }
    *str = '\0';
    return string;
}

void addToPrint(elem_size_t *decArray, uint64_t *size, elem_size_t summand) {
    uint64_t tmp = summand;
    for (uint64_t i = 0; i < *size; i++) {
        if (tmp == 0) {
            return;
        }
        *(decArray + i) = (tmp += *(decArray + i)) % POWER10_9;
        tmp /= POWER10_9;
    }
    if (tmp) {
        *(decArray + *size) = tmp;
        (*size)++;
    }
}

elem_size_t *mul16(elem_size_t *array, uint64_t *size) {
    uint64_t tmp = 0;
    uint64_t max = 16;
    uint64_t size2 = *size;
    for (uint64_t i = 0; i < *size; i++) {
        *(array + i) = (tmp += *(array + i) * max) % POWER10_9;
        tmp /= POWER10_9;
    }
    while (tmp) {
        (*size)++;
        size2++;
        array = realloc(array, (size2) * sizeof(elem_size_t));
        *(array + *size - 1) = tmp % POWER10_9;
        tmp /= POWER10_9;
    }
    return array;
}

uint64_t countBits(elem_size_t *array, uint64_t *size) {
    char *toCheck = malloc(11);
    if (toCheck == 0) {
        fprintf(stderr, "Couldn't allocate memory for a tmp");
        exit(1);
    }
    uint64_t size2 = *size - 1;
    uint16_t addlen = snprintf(toCheck, 11, "%u", *(array + size2));
    uint64_t bits = addlen;
    free(toCheck);
    if (size2 == 0) {
        return bits;
    }
    size2--;
    while (1) {
        if (*(array + size2) == 1000000000) {
            bits += 10;
        } else {
            bits += 9;
        }
        if (size2 != 0) {
            size2--;
        } else {
            break;
        }
    }
    return bits;
}

char *toHex(bignum *a, uint64_t number) {
    uint64_t count = a->size;
    if (count == 0) {
        char *zero = malloc(2);
        if (zero == NULL) {
            fprintf(stderr, "Couldn't allocate memory for a zerostring");
            exit(1);
        }
        *zero = '0';
        *(zero + 1) = 0;
        return zero;
    }
    uint64_t size = count * LOG10_2 * 32 / 9 + 1;
    count--;
    elem_size_t *decArray;
    decArray = malloc(size * sizeof(elem_size_t));
    if (decArray == NULL) {
        fprintf(stderr, "Couldn't allocate memory for decArray in decToPrint");
        exit(1);
    }
    elem_size_t *array = a->array;
    *decArray = 0;
    size = 1;
    uint64_t tmp = *(array + count);
    if (tmp != 0) {
        addToPrint(decArray, &size, tmp);
    }
    if (count > 0) {

        count--;
        while (1) {
            mulToPrint(decArray, &size);
            tmp = *(array + count);
            if (tmp) {
                addToPrint(decArray, &size, tmp);
            }
            if (count == 0) {
                break;
            } else {
                count--;
            }
        }
    }
    char *finalResult, *strf;
    finalResult = strf = malloc(2 * number + 1);
    uint64_t initBits = countBits(decArray, &size);
    uint64_t j = 0;

    uint64_t saveSize = size;
    while (j < number) {

        decArray = mul16(decArray, &size);
        uint64_t newBits = countBits(decArray, &size);
        uint64_t diff;
        if (newBits <= initBits) {
            strf += snprintf(strf, 2, "%x", 0);
            j++;
            continue;
        }
        diff = newBits - initBits;
        char *toCheck = malloc(11);

        if (toCheck == 0) {
            fprintf(stderr, "Couldn't allocate memory for a tmp");
            exit(1);
        }
        uint16_t addlen = snprintf(toCheck, 11, "%u", *(decArray + size - 1));
        free(toCheck);
        if (addlen < diff) {
            // unpleasant case: two digits are in different cells
            elem_size_t tmp;
            tmp = *(decArray + size - 1) * 10 + *(decArray + size - 2) / pow(10, 9);
            strf += snprintf(strf, 2, "%x", tmp);
            *(decArray + size - 2) *= 10;
            *(decArray + size - 2) /= 10;
        } else {
            uint64_t power = 10;
            addlen -= diff;
            power = pow(power, addlen);
            elem_size_t temp = *(decArray + size - 1) / power;
            strf += snprintf(strf, 2, "%x", temp);
            *(decArray + size - 1) %= power;
        }
        size = saveSize;
        j++;
    }
    free(decArray);
    return finalResult;
}

void mulToPrint(elem_size_t *decArray, uint64_t *size) {
    uint64_t tmp = 0;
    uint64_t max = ELEM_SIZE_MAX;
    max++;
    for (uint64_t i = 0; i < *size; i++) {
        *(decArray + i) = (tmp += *(decArray + i) * max) % POWER10_9;
        tmp /= POWER10_9;
    }
    *(decArray + *size) = tmp % POWER10_9;
    (*size)++;
    tmp /= POWER10_9;
    if (tmp) {
        *(decArray + *size) = tmp;
        (*size)++;
    }
}

char *decToPrint(bignum *a) {
    uint64_t count = a->size;
    if (count == 0) {
        char *zero = malloc(2);
        if (zero == NULL) {
            fprintf(stderr, "Couldn't allocate memory for a zerostring\n");
            exit(1);
        }
        *zero = '0';
        *(zero + 1) = 0;
        return zero;
    }
    uint64_t size = count * LOG10_2 * 32 / 9 + 1;
    count--;
    elem_size_t *decArray;
    decArray = malloc(size * sizeof(elem_size_t));
    if (decArray == NULL) {
        fprintf(stderr, "Couldn't allocate memory for decArray in decToPrint\n");
        exit(1);
    }
    elem_size_t *array = a->array;
    *decArray = 0;
    size = 1;
    uint64_t tmp = *(array + count);
    if (tmp != 0) {
        addToPrint(decArray, &size, tmp);
    }
    if (count > 0) {

        count--;
        while (1) {
            mulToPrint(decArray, &size);
            tmp = *(array + count);
            if (tmp) {
                addToPrint(decArray, &size, tmp);
            }
            if (count == 0) {
                break;
            } else {
                count--;
            }
        }
    }
    tmp = *(decArray + size - 1);
    elem_size_t temp = tmp;
    char *string, *str;
    char *toCheck = malloc(11);
    if (toCheck == 0) {
        fprintf(stderr, "Couldn't allocate memory for a tmp\n");
        exit(1);
    }
    // always add plus 1 byte because snprintf always adds and of the string, but we overwrite it then
    uint16_t addLen = snprintf(toCheck, 11, "%u", temp);
    free(toCheck);
    size--;
    uint64_t notOF = size;
    uint64_t max = 10 * notOF + addLen + 1;
    string = str = malloc(max);

    if (string == NULL) {
        fprintf(stderr, "Couldn't allocate memory for string in decToPrint\n");
        exit(1);
    }
    snprintf(str, addLen + 1, "%u", temp);
    str += addLen;
    if (size == 0) {
        *str = '\0';
        free(decArray);
        return string;
    }
    size--;
    while (1) {
        temp = *(decArray + size);
        str += snprintf(str, 11, "%09u", temp);
        if (size == 0) {
            break;
        } else {
            size--;
        }
    }
    *str = '\0';
    free(decArray);
    return string;
}

/**
 *
 * Frees the memory
 *
 */

void freeBigNum(bignum *toFree) {
    if (toFree != NULL) {
        if (toFree->size > 0) {
            free(toFree->array);
        }
        free(toFree);
    }
}

void freeMatrix(matrix *toFree) {
    if (toFree != NULL) {
        freeBigNum(toFree->xn);
        freeBigNum(toFree->xnp1);
        freeBigNum(toFree->xnm1);
        free(toFree);
    }
}

bignum *divideLongDivision(bignum *dividend, bignum *divisor) {
    uint64_t highestBitDividend = calculateHighestBit(dividend->array[dividend->size - 1]);
    uint64_t highestBitDivisor = calculateHighestBit(divisor->array[divisor->size - 1]);
    uint64_t k = N * (dividend->size - 1) + highestBitDividend + 1;
    uint64_t l = N * (divisor->size - 1) + highestBitDivisor + 1;
    bignum *d = new_bignum(divisor->size);
    bignum *r = copy(dividend, l - 1);
    bignum *q = new_bignum(dividend->size);
    bignum *arrayTemp = new_bignum(r->size);
    for (uint64_t i = 0; i <= k - l; i++) {
        uint64_t arrayNumber = dividend->size - 1 - ((31 - highestBitDividend + i + l - 1) / N);
        int64_t bitInArray;
        if (arrayNumber == dividend->size - 1) {

            bitInArray = highestBitDividend - (i + l - 1) % N;
        } else {
            bitInArray = 32 - (i + l - 1) % N;
        }
        uint64_t alpha = dividend->array[arrayNumber] & (1 << bitInArray);
        memcpy(arrayTemp->array, r->array, r->size + 1 * sizeof(elem_size_t));
        arrayTemp->size = r->size;
        arrayTemp = bitShiftLeft(r, 1);
        d->array = arrayTemp->array;
        if (alpha > 0) {
            d = addIntToBignum(d, 1);
        }
        int subtract = compareBignum(d, divisor);
        if (subtract < 1) {
            r = sub(d, divisor);
        } else {
            memcpy(r->array, d->array, (d->size) * sizeof(elem_size_t));
        }
        bitShiftLeft(q, 1);
        if (subtract < 1) {
            q = addIntToBignum(q, 1);
        }
    }
    freeBigNum(d);
    freeBigNum(r);
    return q;
}

bignum *bitShiftRight(bignum *n, uint64_t count) {
    if (n->size == 0 && n->array[0] == 0) {
        return n;
    }
    while (1) {
        for (uint64_t i = 0; i < n->size; i++) {
            n->array[i] = n->array[i] >> 1;
            if (i != n->size - 1 && n->array[i + 1] & 1) {
                n->array[i] |= 1 << (N - 1);
            }
        }
        if (count == 0) {
            break;
        }
        count--;
    }
    zeroJustify(n);
    return n;
}

bignum *bitShiftLeft(bignum *n, uint64_t count) {
    if (n->size == 0 && n->array[0] == 0) {
        return n;
    }
    while (1) {
        if (count == 0) {
            break;
        }
        if (count > 31) {
            arrayShift(n, count / N);
            count = count % N;
        } else {
            uint64_t i = n->size - 1;
            while (1) {
                if (n->array[i] >= TWO_POW_31) {
                    if (i == n->size - 1) {
                        n->array = realloc(n->array, (n->size + 1) * sizeof(elem_size_t));
                        n->size++;
                        n->array[n->size - 1] = 0;
                    }
                    n->array[i + 1] += 1;
                }
                n->array[i] = n->array[i] << 1;
                if (i == 0) {
                    break;
                }
                i--;
            }
            count--;
        }
    }
    return n;
}

bignum *addIntToBignum(bignum *n, uint64_t toAdd) {
    bignum *temp = new_bignum(1);
    temp->array[0] = toAdd;
    return add(n, temp);
}

bignum *copy(bignum *from, uint64_t countBits) {
    uint64_t size = (countBits - 1) / N;
    bignum *res = new_bignum(size + 1);
    uint64_t currentArrayFrom = from->size - 1;
    uint64_t currentBitFrom = calculateHighestBit(from->array[from->size - 1]);
    uint64_t currentArraySet = (countBits - 1) / N;
    uint64_t currentBitSet = (countBits - 1) % N;
    uint64_t counter = countBits - 1;
    while (1) {
        if (from->array[currentArrayFrom] & 1ULL << currentBitFrom) {
            res->array[currentArraySet] |= 1ULL << currentBitSet;
        }
        if (currentBitFrom == 0) {
            currentArrayFrom--;
            currentBitFrom = 31;
        } else {
            currentBitFrom--;
        }
        if (currentBitSet == 0) {
            currentArraySet--;
            currentBitSet = 31;
        } else {
            currentBitSet--;
        }
        if (counter == 0) {
            break;
        }
        counter--;
    }
    return res;
}


