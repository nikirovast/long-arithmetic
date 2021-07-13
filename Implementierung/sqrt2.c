#include "assert.h"
#include "errno.h"
#include "math.h"
#include "time.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#pragma warning(disable : 4996)

#define N 32
#define TWO_POW_31 2147483648
#define ELEM_SIZE_MAX UINT32_MAX
#define LOG2_10 3.3
#define LOG10_2 0.3
#define POWER10_9 1000000000
typedef uint32_t elem_size_t;
typedef struct bigNum bigNum;
typedef struct matrix matrix;
static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};
struct bigNum
{
    size_t size;
    elem_size_t *array;
};

extern void sum(uint64_t n, bigNum *xn, bigNum *xnp1);
bigNum *mul(bigNum *xn, bigNum *xnp1);
void zero_justify(bigNum *n);
void array_shift(bigNum *n, int count);
int compare_bignum(bigNum *xn, bigNum *xnp1);
bigNum *div2bignums(bigNum *xn, bigNum *xnp1);
bigNum *add(bigNum *xn, bigNum *xnp1);
bigNum *sub(bigNum *xn, bigNum *xnp1);
bigNum *mul2num(elem_size_t x, elem_size_t y);
matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2);
matrix *matrixBinaryExponentiation(unsigned long long n, int highestBit);
int calculateHighestBit(unsigned long long n);
uint64_t convertAccToN(uint64_t numDigits);
char *hexToPrint(bigNum *a);
char *decToPrint(bigNum *a);
void printUsage(char **argv);
bigNum *divideLongDivision(bigNum *dividend, bigNum *divisor);
bigNum *bitShiftRight(bigNum *n, int count);
bigNum *bitShiftLeft(bigNum *n, int count);
bigNum *addIntToBigNum(bigNum *n, int toAdd);

bigNum *new_bignum(size_t size)
{
    bigNum *res = malloc(sizeof(bigNum));
    if (res == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for a struct");
        exit(1);
    }
    res->size = size;
    res->array = malloc(sizeof(elem_size_t) * size);
    if (res->array == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for an array");
        exit(1);
    }
    return res;
}

struct matrix
{
    bigNum *xn;
    bigNum *xnp1;
    bigNum *xnm1;
};

// in VSCode 1st is in RCX, 2nd in RDX, 3rd in R8
int main(int argc, char **argv) {
    // printf("%i", convertAccToN(10));
    // bigNum* a = new_bignum(2);
    // bigNum* b = new_bignum(1);
    //*(a->array + 0) = 4294967295;
    //*(a->array + 1) = 1;
    //*(b->array + 0) = 4294967295;
    // char* r = hexPrint(a);
    // printf("%u, %u, %u", a->array[0], a->array[1], b->array[0]);
    // bigNum* res = mul(a, b);

    // uint64_t s = res->size;

    // printf("\n");
    // printf("Size: %lu\n", res->size);
    // for (int i = res->size - 1; i >= 0; i--) {
    //  printf("%u\n", res->array[i]);
    //}

    //// test mul1
    // bigNum* first = new_bignum(2);
    // bigNum* second = new_bignum(1);
    //*(first->array + 0) = 1234;
    //*(first->array + 1) = 1;
    //*(second->array + 0) = 129;
    // bigNum* res1 = mul(first, second);
    //// 554050940370
    // for (int i = res1->size - 1; i >= 0; i--) {
    //  printf("%u\n", res1->array[i]);
    //}
    // assert(res1->array[0] == 159186);
    // assert(res1->array[1] == 129);

    //// test mul2
    // bigNum* third = new_bignum(1);
    // bigNum* fourth = new_bignum(1);
    //*(third->array) = 65537;
    //*(fourth->array) = 89340;
    // bigNum* res2 = mul(third, fourth);
    //// 5855075580
    // for (int i = res2->size - 1; i >= 0; i--) {
    //  printf("%u\n", res2->array[i]);
    //}
    // assert(res2->array[0] == 1560108284);
    // assert(res2->array[1] == 1);

    //// test div1
    // bigNum* divisor = new_bignum(2);
    // bigNum* dividend = new_bignum(1);
    //*(divisor->array) = 90000;
    //*(divisor->array + 1) = 1;
    //*(dividend->array) = 7500;
    // bigNum* res3 = div2bignums(divisor, dividend);
    // printf("Result of division: \n");
    // printf("Size: %lu\n", res3->size);
    // for (long i = res3->size - 1; i >= 0; i--) {
    //  printf("%u\n", res3->array[i]);
    //}

    /*struct matrix *matrix = malloc(3 * sizeof(bigNum));
    if (matrix == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for a matrix");
        exit(1);
    }
    matrix->xnp1 = new_bignum(1);
    matrix->xnp1->array[0] = 2;
    matrix->xn = new_bignum(1);
    matrix->xn->array[0] = 1;
    matrix->xnm1 = new_bignum(1);
    matrix->xnm1->array[0] = 0;

    int n = 20;
    struct matrix *res4;
    struct matrix *res5;
    int op = convertAccToN(n);
    res5 = matrixBinaryExponentiation(matrix, op, calculateHighestBit(op));
     */
    // array_shift(res5->xn, 1);
    // res4 = div2bignums(res5->xn, res5->xnp1);
    //char *result = hexToPrint(res5->xn);
    //char *test2 = decToPrint(res5->xn);
    // double result = (double)res5->xn->array[0] / (double)res5->xnp1->array[0];

    // char *result = hexToPrint(res4->xn);

    // bigNum *a = new_bignum(2);
    // bigNum *b = new_bignum(2);
    //*(a->array) = 3521874100;
    //*(a->array + 1) = 2154513721;
    ////*a->array + 2) = 1774;
    //*(b->array) = 3521874100;
    //*(b->array + 1) = 2154513721;
    ////*(b->array + 2) = 1774;
    // bigNum *res = mul(a, b);
    // char *result = hexToPrint(res);
    // printf("End matrix: %u, %u, %u with size of xnp1 %lu",
    // res4->xnm1->array[0],
    //       res4->xn->array[0], res4->xnp1->array[0], res4->xnp1->size);
    int hexadecimal = 0;
    int c;
    if (argc == 1 || argc > 3) {
        printUsage(argv);
        return 1;
    }
    while (1) {
        int option_index = 0;
        c = getopt_long(argc, argv, "hx", long_options, &option_index);
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
            default: {
                fprintf(stderr, "Unknown option was used");
                exit(1);
            }
        }
    }
    // well known bug (or feature) with getopt_long changing the order of arguments for some reason
    uint64_t n = strtoull(argv[1 + hexadecimal], NULL, 0);
    if (errno == ERANGE) {
        fprintf(stderr, "the given number can not be represented, please pick a number < UINT64_MAX");
        return 1;
    }
    if (n == 0ULL || *argv[1 + hexadecimal] == '-') {
        fprintf(stderr, "invalid number of iterations: it should be > 0 or the given parameter was not a number");
        return 1;
    }
    struct matrix *res5;
    uint64_t op = convertAccToN(n);
    bigNum *help = new_bignum(1);
    /*help->array[0] = 0b1100011000100000000000000000000;
    help->array[1] = 0b1101011110001110101111000101101;
    help->array[2] = 0b101;
     */
    help->array[0] = 10;
    res5 = matrixBinaryExponentiation(op, calculateHighestBit(op));

    if (n < 1) {
        double result = (double) res5->xn->array[0] / (double) res5->xnp1->array[0];
        printf("Result of division after %llu operations: %.15f", op, result);
    } else {


        bigNum *helper = mul(res5->xn, help);
        //bigNum *res = div2bignums(helper, res5->xnp1);
        /*bigNum* first = new_bignum(1);
        bigNum* second = new_bignum(1);
        first->array[0] = 185;
        second->array[0] = 13;
        bigNum *res = divideLongDivision(first, second);
         */
        bigNum *res = divideLongDivision(helper, res5->xnp1);
        printf("Result of division after %llu operations: %u\n", op , res->array[0]);
        printf("%u\n", res->array[1]);
        printf("%u", res->array[2]);
        return 0;
    }
}

/**
 *
 * Multiply two matrices.
 *
 */

void printUsage(char **argv)
{
    printf("Usage: %s <number of iterations> \n", argv[0]);
    printf("-x: sets output numeral system to hexadecimal, default: decimal\n");
    printf("-h|--help: for usage help");
}

matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2)
{
    matrix *res = malloc(2 * sizeof(*matrix2));
    if (res == 0)
    {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixMul");
        exit(1);
    }
    res->xnm1 = add(mul(matrix1->xnm1, matrix2->xnm1), mul(matrix1->xn, matrix2->xn));
    res->xn = add(mul(matrix1->xnm1, matrix2->xn), mul(matrix1->xn, matrix2->xnp1));
    res->xnp1 = add(mul(matrix1->xn, matrix2->xn), mul(matrix1->xnp1, matrix2->xnp1));
    return res;
}

/**
 *
 * Binary exponentiation of the matrix to the exponent n as described in pdf
 *
 */

matrix *matrixBinaryExponentiation(unsigned long long n, int highestBit)
{
    struct matrix *matrix = malloc(3 * sizeof (struct bigNum));
    matrix->xn = new_bignum(1);
    matrix->xn->array[0] = 1;
    matrix->xnm1 = new_bignum(1);
    matrix->xnm1->array[0] = 0;
    matrix->xnp1 = new_bignum(1);
    matrix->xnp1->array[0] = 2;
    struct matrix *matrixInitial = malloc(3 * sizeof(bigNum));
    if (matrixInitial == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixBinaryExp");
        exit(1);
    }
    matrixInitial->xnm1 = new_bignum(1);
    *(matrixInitial->xnm1->array) = 0;
    matrixInitial->xn = new_bignum(1);
    *(matrixInitial->xn->array) = 1;
    matrixInitial->xnp1 = new_bignum(1);
    *(matrixInitial->xnp1->array) = 2;
    for (int i = highestBit - 1; i >= 0; i--)
    {
        matrix = matrixMultiplication(matrix, matrix);
        if (n >> i & 1)
        {
            matrix = matrixMultiplication(matrix, matrixInitial);
        }
    }
    free(matrixInitial->xn);
    free(matrixInitial->xnm1);
    free(matrixInitial->xnp1);
    free(matrixInitial);
    return matrix;
}

/**
 *
 * Calculate highest bit of the number
 *
 */

int calculateHighestBit(unsigned long long n)
{
    int order = -1;
    for (int i = 0; i < N; i++)
    {
        if ((n >> i) & 1)
        {
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

bigNum *add(bigNum *xn, bigNum *xnp1)
{
    int flag = 0;
    if (xn->size > xnp1->size)
    {
        bigNum *tmp = xn;
        xn = xnp1;
        xnp1 = tmp;
        flag = 1;
    }

    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;
    bigNum *res = new_bignum(size2);
    int i;
    elem_size_t carry = 0;
    for (i = 0; i < size1; i++)
    {
        elem_size_t a = *(xnp1->array + i);
        elem_size_t b = *(xn->array + i);
        *(res->array + i) = a + b + carry;
        carry = (a > ELEM_SIZE_MAX - b);
    }
    while (i < size2)
    {
        elem_size_t a = *(xnp1->array + i);
        *(res->array + i) = a + carry;
        carry = (a > ELEM_SIZE_MAX - carry);
        i++;
    }
    if (carry)
    {
        elem_size_t *tmp = realloc(res->array, sizeof(elem_size_t) * (size2 + 1));
        if (tmp == NULL)
        {
            fprintf(stderr, "Couldn't reallocate memory for a carry in sum");
            exit(1);
        }
        else
        {
            res->array = tmp;
            res->size += 1;
            *(res->array + i) = carry;
        }
    }
    if (flag != 0)
    {
        bigNum *tmp = xn;
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

bigNum *sub(bigNum *xn, bigNum *xnp1)
{
    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;
    bigNum *res = new_bignum(sizeof(elem_size_t) * size1);
    res->size = size1;
    int i;
    elem_size_t carry = 0;
    for (i = 0; i < size2; i++)
    {
        elem_size_t a = *(xn->array + i);
        elem_size_t b = *(xnp1->array + i);
        if ((a < 1 && (b > 0 || carry > 0)) || a - carry < b)
        {
            *(res->array + i) = ELEM_SIZE_MAX - b - carry + a + 1;
            carry = 1;
        }
        else
        {
            *(res->array + i) = a - b - carry;
            carry = 0;
        }
    }
    while (i < size1)
    {
        elem_size_t a = *(xn->array + i);
        if (a >= carry)
        {
            *(res->array + i) = a - carry;
            carry = 0;
        }
        else
        {
            *(res->array + i) = ELEM_SIZE_MAX - carry + a + 1;
            carry = 1;
        }
        i++;
    }
    zero_justify(res);
    return res;
}

/**
 *
 * Multiply two big numbers xnp1 by xn.
 *
 */

bigNum *mul(bigNum *xn, bigNum *xnp1)
{
    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;

    if (xn->size == 0 || xnp1->size == 0)
    {
        bigNum *res = malloc(sizeof(*xn));
        if (res == 0)
        {
            fprintf(stderr, "Couldn't allocate memory for a zero res in mul");
            exit(1);
        }
        res->size = 0;
        return res;
    }
    else if (size1 > 1 || size2 > 1)
    {
        uint64_t aNewSize = size1 / 2;
        uint64_t bNewSize = size1 - aNewSize;
        bigNum *b = new_bignum(bNewSize);
        if (b == NULL)
        {
            fprintf(stderr, "Couldn't allocate memory for b in mul");
            exit(1);
        }
        b->size = bNewSize;
        // we need to write in b the second part of xn
        // void *memcpy(void *dest, const void *src, std::size_t count);
        memcpy(b->array, xn->array, bNewSize * sizeof(elem_size_t));
        // b->array = xn->array;

        bigNum *a = new_bignum(aNewSize);

        if (a == NULL)
        {
            fprintf(stderr, "Couldn't allocate memory for a in mul");
            exit(1);
        }
        a->size = aNewSize;
        // memcpy(a->array, (xn->array + bNewSize), aNewSize * sizeof(elem_size_t));
        a->array = xn->array + bNewSize;

        uint64_t cNewSize = size2 / 2;
        uint64_t dNewSize = size2 - cNewSize;
        bigNum *c = new_bignum(cNewSize);
        if (c == NULL)
        {
            fprintf(stderr, "Couldn't allocate memory for c in mul");
            exit(1);
        }
        c->size = cNewSize;
        c->array = xnp1->array + dNewSize;
        // memcpy(c->array, (xnp1->array + dNewSize), cNewSize * sizeof(elem_size_t));

        bigNum *d = new_bignum(dNewSize);
        if (d == NULL)
        {
            fprintf(stderr, "Couldn't allocate memory for d in mul");
            exit(1);
        }
        d->size = dNewSize;
        memcpy(d->array, xnp1->array, dNewSize * sizeof(elem_size_t));
        // d->array = xnp1->array;
        // memcpy(d->array, xnp1->array, xnp1->size * sizeof(unsigned long long));
        bigNum *ac = mul(a, c);
        bigNum *bd = mul(b, d);
        bigNum *a_add_b = add(a, b);
        bigNum *c_add_d = add(c, d);

        bigNum *abcd = mul(a_add_b, c_add_d);
        free(a_add_b);
        free(c_add_d);
        bigNum *abcd_sub_ac = sub(abcd, ac);
        bigNum *ad_add_bc = sub(abcd_sub_ac, bd);
        free(a);
        free(b);
        free(c);
        free(d);
        if (ac->size != 0)
        {
            elem_size_t *tmp = realloc(ac->array, sizeof(elem_size_t) * (ac->size + 2));
            if (tmp == NULL)
            {
                fprintf(stderr, "Couldn't reallocate memory for a shift in sum");
                exit(1);
            }
            else
            {
                ac->array = tmp;
                ac->size += 2;
                for (int i = ac->size - 1; i > 1; i--)
                {
                    *(ac->array + i) = *(ac->array + i - 2);
                }
                *(ac->array) = 0;
                *(ac->array + 1) = 0;
            }
        }
        if (ad_add_bc->size != 0)
        {
            elem_size_t *tmp = realloc(ad_add_bc->array, sizeof(elem_size_t) * (ad_add_bc->size + 1));
            if (tmp == NULL)
            {
                fprintf(stderr, "Couldn't allocate memory for a zero shift in sum");
                exit(1);
            }
            else
            {
                ad_add_bc->array = tmp;
                ad_add_bc->size += 1;

                for (int i = ad_add_bc->size - 1; i > 0; i--)
                {
                    *(ad_add_bc->array + i) = *(ad_add_bc->array + i - 1);
                }
                *(ad_add_bc->array) = 0;
            }
        }
        bigNum *res = add(add(ac, ad_add_bc), bd);
        free(abcd);
        free(abcd_sub_ac);
        free(ad_add_bc);
        zero_justify(res);

        return res;
    }
    else
    {
        return mul2num(*(xn->array), *(xnp1->array));
    }
}

/**
 *
 * Multiply two elem_size_t numbers x by y.
 *
 */

bigNum *mul2num(elem_size_t x, elem_size_t y)
{
    if (x == 0 || y == 0)
    {
        bigNum *res = malloc(sizeof(*res));
        if (res == 0)
        {
            fprintf(stderr, "Couldn't allocate memory for a zero res in mul");
            exit(1);
        }
        res->size = 0;
        return res;
    }
    else
    {
        size_t shift_size = N / 2;
        elem_size_t a = x >> shift_size;
        elem_size_t b = (x << shift_size) >> shift_size;
        elem_size_t c = y >> shift_size;
        elem_size_t d = (y << shift_size) >> shift_size;
        bigNum ac;
        ac.size = 2;
        elem_size_t ac_array[] = {0, a * c};

        ac.array = ac_array;

        bigNum ad;
        ad.size = 2;
        elem_size_t ad_array[] = {((a * d) << shift_size), (a * d) >> shift_size};
        ad.array = ad_array;

        bigNum bc;
        bc.size = 2;
        elem_size_t bc_array[] = {(b * c) << shift_size, (b * c) >> shift_size};
        bc.array = bc_array;

        bigNum bd;
        bd.size = 1;
        elem_size_t bd_array[] = {b * d};
        bd.array = bd_array;
        bigNum *ac_add_ad = add(&ac, &ad);
        bigNum *ac_add_ad_add_bc = add(ac_add_ad, &bc);
        free(ac_add_ad);
        bigNum *res = add(ac_add_ad_add_bc, &bd);
        zero_justify(res);
        return res;
    }
}

/**
 *
 * Divide big number xn by big number xnp1.
 *
 */

bigNum *div2bignums(bigNum *xn, bigNum *xnp1)
{
    size_t size = xn->size;
    bigNum *res = new_bignum(size);
    bigNum *row = new_bignum(1);
    uint64_t t = xnp1->array[xnp1->size - 1] * 2;
    for (long i = xn->size - 1; i >= 0; i--)
    {
        array_shift(row, 1);
        row->array[0] = xn->array[i];
        res->array[i] = 0;
        while (compare_bignum(row, xnp1) < 1) {
            if (t < row->array[xnp1->size-1]) {
                bigNum *temp = new_bignum(1);
                temp->array[0] = row->array[xnp1->size-1] / t;
                row = sub(row, mul(xnp1, temp));
                res->array[i] += temp->array[0];

            }
            else {
                res->array[i]++;
                row = sub(row, xnp1);
            }
            zero_justify(row);
        }
    }
    zero_justify(res);
    free(row);
    return res;
}

int compare_bignum(bigNum *xn, bigNum *xnp1)
{
    if (xn->size > xnp1->size)
    {
        return -1;
    }
    else if (xnp1->size > xn->size)
    {
        return 1;
    }
    for (long i = xnp1->size - 1; i >= 0; i--)
    {
        if (xn->array[i] > xnp1->array[i])
        {
            return -1;
        }
        if (xn->array[i] < xnp1->array[i])
        {
            return 1;
        }
    }
    return 0;
}

void array_shift(bigNum *n, int count)
{
    long i;
    elem_size_t *tmp = realloc(n->array, (count + n->size) * sizeof(elem_size_t));
    if (tmp == NULL)
    {
        fprintf(stderr, "Couldn't reallocate memory for a n in arrayShift");
        exit(1);
    }
    else
    {
        n->array = tmp;
    }
    if ((n->size == 1) && (n->array[0] == 0))
        return;

    for (i = n->size - 1; i >= 0; i--) {
        n->array[i + count] = n->array[i];
    }

    for (i = 0; i < count; i++) {
        *(n->array + i) = 0;
    }
    n->size = n->size + count;
}

void zero_justify(bigNum *n)
{
    while ((n->size > 1) && (n->array[n->size - 1] == 0))
    {
        n->size--;
    }
}

/**
 *
 * Counts num of iteratons needed to achieve certain accuracy
 *
 */

uint64_t convertAccToN(uint64_t numDigits)
{
    return (uint64_t)ceil(1 + (LOG2_10 / 2) * numDigits);
}

/**
 *
 * Converts bigNum to hexadecimal format string
 *
 */

char *hexToPrint(bigNum *a)
{
    uint64_t count = a->size;
    int numCharinOneElem = N / 4;
    // One extra for a string terminator
    int len = count * numCharinOneElem + 1;
    count--;
    elem_size_t *array = a->array;

    elem_size_t first = *(array + count);
    if (first <= 0xf)
    {
        len -= 7;
    }
    else if (first <= 0xff)
    {
        len -= 6;
    }
    else if (first <= 0xfff)
    {
        len -= 5;
    }
    else if (first <= 0xffff)
    {
        len -= 4;
    }
    else if (first <= 0xfffff)
    {
        len -= 3;
    }
    else if (first <= 0xffffff)
    {
        len -= 2;
    }
    else if (first <= 0xfffffff)
    {
        len--;
    }

    char *string, *str;
    str = string = malloc(len);
    if (string == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for a string in hexPrint");
        exit(1);
    }

    int k = sprintf(str, "%x", first);
    str = str + k;
    count--;
    for (; count >= 0; count--)
    {
        elem_size_t c = *(array + count);
        if (c <= 0xf)
        {
            sprintf(str, "0000000%x", c);
        }
        else if (c <= 0xff)
        {
            sprintf(str, "000000%x", c);
        }
        else if (c <= 0xfff)
        {
            sprintf(str, "00000%x", c);
        }
        else if (c <= 0xffff)
        {
            sprintf(str, "0000%x", c);
        }
        else if (c <= 0xfffff)
        {
            sprintf(str, "000%x", c);
        }
        else if (c <= 0xffffff)
        {
            sprintf(str, "00%x", c);
        }
        else if (c <= 0xfffffff)
        {
            sprintf(str, "0%x", c);
        }
        else
        {
            sprintf(str, "%x", c);
        }
        str += numCharinOneElem;
    }
    *str = '\0';
    return string;
}

void addToPrint(elem_size_t *decArray, int *size, elem_size_t summand)
{
    uint64_t tmp = summand;
    for (int i = 0; i < *size; i++)
    {
        if (!tmp)
            return;
        *(decArray + i) = (tmp += *(decArray + i)) % POWER10_9;
        tmp /= POWER10_9;
    }
    if (tmp)
    {
        *(decArray + *size) = tmp;
        (*size)++;
    }
}

void mulToPrint(elem_size_t *decArray, int *size)
{
    uint64_t tmp = 0;
    for (int i = 0; i < *size; i++)
    {
        *(decArray + i) = (tmp += *(decArray + i) * (ELEM_SIZE_MAX + 1)) % POWER10_9;
        tmp /= POWER10_9;
    }
    *(decArray + *size) = tmp % POWER10_9;
    (*size)++;
    tmp /= POWER10_9;
    if (tmp)
    {
        *(decArray + *size) = tmp;
        (*size)++;
    }
}

char *decToPrint(bigNum *a)
{
    int count = a->size;
    int size = count * LOG10_2 * 32 / 9 + 1;
    count--;
    elem_size_t *decArray;
    decArray = malloc(size * sizeof(elem_size_t));
    if (decArray == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for decArray in decToPrint");
        exit(1);
    }
    elem_size_t *bytes = a->array;
    *decArray = 0;
    size = 1;
    uint64_t tmp = *(bytes + count);
    if (tmp)
    {
        addToPrint(decArray, &size, tmp);
    }

    while (count--)
    {
        mulToPrint(decArray, &size);
        tmp = *(bytes + count);
        if (tmp)
            addToPrint(decArray, &size, tmp);
    }

    tmp = *(decArray + size - 1);
    char *string, *str;
    int len = 9;
    if (tmp < 10)
    {
        len = 1;
    }
    else if (tmp < 100)
    {
        len = 2;
    }
    else if (tmp < 1000)
    {
        len = 3;
    }
    else if (tmp < 10000)
    {
        len = 4;
    }
    else if (tmp < 100000)
    {
        len = 5;
    }
    else if (tmp < 1000000)
    {
        len = 6;
    }
    else if (tmp < 10000000)
    {
        len = 7;
    }
    else if (tmp < 100000000)
    {
        len = 8;
    }
    string = str = malloc(9 * size + len + 1);

    if (string == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for string in decToPrint");
        exit(1);
    }
    str += sprintf(str, "%llu", tmp);
    size--;
    while (size--)
    {
        tmp = *(decArray + size);
        if (tmp < 10)
        {
            sprintf(str, "00000000%llu", tmp);
        }
        else if (tmp < 100)
        {
            sprintf(str, "0000000%llu", tmp);
        }
        else if (tmp < 1000)
        {
            sprintf(str, "000000%llu", tmp);
        }
        else if (tmp < 10000)
        {
            sprintf(str, "00000%llu", tmp);
        }
        else if (tmp < 100000)
        {
            sprintf(str, "0000%llu", tmp);
        }
        else if (tmp < 1000000)
        {
            sprintf(str, "000%llu", tmp);
        }
        else if (tmp < 10000000)
        {
            sprintf(str, "00%llu", tmp);
        }
        else if (tmp < 100000000)
        {
            sprintf(str, "0%llu", tmp);
        }
        else
        {
            sprintf(str, "%llu", tmp);
        }

        str += 9;
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

void freeBigNum(bigNum *toFree)
{
    if (toFree->array != NULL)
    {
        free(toFree->array);
    }
    free(toFree);
}

bigNum *divideLongDivision(bigNum *dividend, bigNum *divisor) {
    int highestBitDividend = calculateHighestBit(dividend->array[dividend->size - 1]);
    int highestBitDivisor = calculateHighestBit(divisor->array[divisor->size - 1]);
    uint64_t k = N * (dividend->size - 1) + highestBitDividend + 1;
    uint64_t l = N * (divisor->size - 1) + highestBitDivisor + 1;
    bigNum *r = new_bignum(divisor->size);
    bigNum *d = new_bignum(divisor->size);
    memcpy(r->array, divisor->array, sizeof(elem_size_t) * (divisor->size));
    //clearing a bit according to long division algorithm
    r->array[r->size-1] &= ~(1UL << highestBitDivisor);
    bigNum *q = new_bignum(dividend->size);
    bigNum *arrayTemp = new_bignum(r->size);
    for(int i = 0; i <= k - l; i++) {
        uint64_t arrayNumber = dividend->size - 1 - ((31 - highestBitDividend + i + l - 1) / N);
        uint64_t bitInArray;
        if (arrayNumber == dividend->size - 1) {
            bitInArray = highestBitDividend - (i + l - 1) % N;
        }
        else {
            bitInArray = 31 - (i + l - 1) % N;
        }
        int alpha = dividend->array[arrayNumber] & (1<<bitInArray);
        memcpy(arrayTemp->array, r->array, (r->size) * sizeof (elem_size_t));
        arrayTemp = bitShiftLeft(r, 1);
        d->array = arrayTemp->array;
        if (alpha > 0) {
            d = addIntToBigNum(d, 1);
        }
        int subtract = compare_bignum(d, divisor);
        if (subtract < 1) {
            r = sub(d, divisor);
        }
        else {
            memcpy(r->array, d->array, (d->size) * sizeof (elem_size_t));
        }
        bitShiftLeft(q, 1);
        if (subtract < 1) {
            q = addIntToBigNum(q, 1);
        }
    }
    free(arrayTemp);
    printf("Quotient: %u and rest: %u\n", q->array[0], r->array[0]);
    return q;
}

bigNum *bitShiftRight(bigNum *n, int count) {
    if (n->size == 1 && n->array[0] == 0) {
        return n;
    }
   while (count > 0) {
       for (int i = 0; i < n->size; i++) {
           n->array[i] = n->array[i]>>1;
           if (i != n->size - 1 && n->array[i+1]&1) {
               n->array[i] |= 1<<(N-1);
           }
       }
       count--;
   }
   zero_justify(n);
   return n;
}

bigNum *bitShiftLeft(bigNum *n, int count) {
    if (n->size == 1 && n->array[0] == 0) {
        return n;
    }
    while (count > 0) {
        for (int i = n->size - 1; i >= 0; i--) {
            if (n->array[i] >= TWO_POW_31) {
                if (i == n->size - 1) {
                    n->array = realloc(n->array, (n->size + 1) * sizeof(elem_size_t));
                    n->size++;
                    n->array[n->size - 1] = 0;
                }
                n->array[i + 1] += 1;
            }
                n->array[i] = n->array[i]<<1;
        }
        count--;
    }
    return n;
}

bigNum *addIntToBigNum(bigNum *n, int toAdd) {
    bigNum *temp = new_bignum(1);
    temp->array[0] = toAdd;
    return add(n, temp);
}
