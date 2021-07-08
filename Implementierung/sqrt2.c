#include "assert.h"
#include "errno.h"
#include "math.h"
#include "time.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma warning(disable : 4996)

#define N 32
#define ELEM_SIZE_MAX UINT32_MAX

typedef uint32_t elem_size_t;
typedef struct bignum bignum;
typedef struct matrix matrix;
struct bignum
{
    size_t size;
    elem_size_t *array;
};

extern void sum(uint64_t n, bignum *xn, bignum *xnp1);
bignum *mul(bignum *xn, bignum *xnp1);
void zero_justify(bignum *n);
void array_shift(bignum *n, int count);
int compare_bignum(bignum *xn, bignum *xnp1);
bignum *div2bignums(bignum *xn, bignum *xnp1);
// void sqrt2(uint64_t n, bignum *xn, bignum *xnp1);#
bignum *add(bignum *xn, bignum *xnp1);
bignum *sub(bignum *xn, bignum *xnp1);
bignum *mul2num(elem_size_t x, elem_size_t y);
matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2);
matrix *matrixBinaryExponentiation(matrix *matrix, unsigned long long n, int highestBit);
int calculateHighestBit(unsigned long long n);
uint64_t convertAccToN(uint64_t numDigits);
char *hexToPrint(bignum *a);

bignum *new_bignum(size_t size)
{
    bignum *res = malloc(sizeof(bignum));
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
    bignum *xn;
    bignum *xnp1;
    bignum *xnm1;
};

// in VSCode 1st is in RCX, 2nd in RDX, 3rd in R8
int main(int argc, char **argv)
{
    // printf("%i", convertAccToN(10));
    // bignum* a = new_bignum(2);
    // bignum* b = new_bignum(1);
    //*(a->array + 0) = 4294967295;
    //*(a->array + 1) = 1;
    //*(b->array + 0) = 4294967295;
    // char* r = hexPrint(a);
    // printf("%u, %u, %u", a->array[0], a->array[1], b->array[0]);
    // bignum* res = mul(a, b);

    // uint64_t s = res->size;

    // printf("\n");
    // printf("Size: %lu\n", res->size);
    // for (int i = res->size - 1; i >= 0; i--) {
    //  printf("%u\n", res->array[i]);
    //}

    //// test mul1
    // bignum* first = new_bignum(2);
    // bignum* second = new_bignum(1);
    //*(first->array + 0) = 1234;
    //*(first->array + 1) = 1;
    //*(second->array + 0) = 129;
    // bignum* res1 = mul(first, second);
    //// 554050940370
    // for (int i = res1->size - 1; i >= 0; i--) {
    //  printf("%u\n", res1->array[i]);
    //}
    // assert(res1->array[0] == 159186);
    // assert(res1->array[1] == 129);

    //// test mul2
    // bignum* third = new_bignum(1);
    // bignum* fourth = new_bignum(1);
    //*(third->array) = 65537;
    //*(fourth->array) = 89340;
    // bignum* res2 = mul(third, fourth);
    //// 5855075580
    // for (int i = res2->size - 1; i >= 0; i--) {
    //  printf("%u\n", res2->array[i]);
    //}
    // assert(res2->array[0] == 1560108284);
    // assert(res2->array[1] == 1);

    //// test div1
    // bignum* divisor = new_bignum(2);
    // bignum* dividend = new_bignum(1);
    //*(divisor->array) = 90000;
    //*(divisor->array + 1) = 1;
    //*(dividend->array) = 7500;
    // bignum* res3 = div2bignums(divisor, dividend);
    // printf("Result of division: \n");
    // printf("Size: %lu\n", res3->size);
    // for (long i = res3->size - 1; i >= 0; i--) {
    //  printf("%u\n", res3->array[i]);
    //}

    struct matrix *matrix = malloc(3 * sizeof(bignum));
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

    int n = 100;
    struct matrix *res4;
    int hb = calculateHighestBit(n);
    res4 = matrixBinaryExponentiation(matrix, n, hb);
    char *result = hexToPrint(res4->xn);

    // bignum *a = new_bignum(2);
    // bignum *b = new_bignum(2);
    //*(a->array) = 3521874100;
    //*(a->array + 1) = 2154513721;
    ////*a->array + 2) = 1774;
    //*(b->array) = 3521874100;
    //*(b->array + 1) = 2154513721;
    ////*(b->array + 2) = 1774;
    // bignum *res = mul(a, b);
    // char *result = hexToPrint(res);
    // printf("End matrix: %u, %u, %u with size of xnp1 %lu",
    // res4->xnm1->array[0],
    //       res4->xn->array[0], res4->xnp1->array[0], res4->xnp1->size);
    return 0;
}

/**
 *
 * Multiply two matrices.
 *
 */

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
    // useless
    // zero_justify(res->xnm1);
    // zero_justify(res->xn);
    // zero_justify(res->xnp1);
    return res;
}

/**
 *
 * Binary exponentiation of the matrix to the exponent n as described in pdf
 *
 */

matrix *matrixBinaryExponentiation(matrix *matrix, unsigned long long n, int highestBit)
{
    struct matrix *matrixInitial = malloc(3 * sizeof(bignum));
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

bignum *add(bignum *xn, bignum *xnp1)
{
    int flag = 0;
    if (xn->size > xnp1->size)
    {
        bignum *tmp = xn;
        xn = xnp1;
        xnp1 = tmp;
        flag = 1;
    }

    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;
    bignum *res = new_bignum(sizeof(elem_size_t) * size2);
    res->size = size2;
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

bignum *sub(bignum *xn, bignum *xnp1)
{
    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;
    bignum *res = new_bignum(sizeof(elem_size_t) * size1);
    res->size = size1;
    int i;
    elem_size_t carry = 0;
    for (i = 0; i < size2; i++)
    {
        elem_size_t a = *(xn->array + i);
        elem_size_t b = *(xnp1->array + i);
        if (a < 1 && (b > 0 || carry > 0) || a - carry < b)
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

bignum *mul(bignum *xn, bignum *xnp1)
{
    uint64_t size1 = xn->size;
    uint64_t size2 = xnp1->size;

    if (xn->size == 0 || xnp1->size == 0)
    {
        bignum *res = malloc(sizeof(*xn));
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
        bignum *b = new_bignum(bNewSize);
        b->size = bNewSize;
        // we need to write in b the second part of xn
        // void *memcpy(void *dest, const void *src, std::size_t count);
        memcpy(b->array, xn->array, bNewSize * sizeof(elem_size_t));
        // b->array = xn->array;

        bignum *a = new_bignum(aNewSize);
        a->size = aNewSize;
        memcpy(a->array, (xn->array + bNewSize), aNewSize * sizeof(elem_size_t));
        // a->array = xn->array + bNewSize;

        uint64_t cNewSize = size2 / 2;
        uint64_t dNewSize = size2 - cNewSize;
        bignum *c = new_bignum(cNewSize);
        c->size = cNewSize;
        // c->array = xnp1->array + dNewSize;
        memcpy(c->array, (xnp1->array + dNewSize), cNewSize * sizeof(elem_size_t));

        bignum *d = new_bignum(dNewSize);
        d->size = dNewSize;
        memcpy(d->array, xnp1->array, dNewSize * sizeof(elem_size_t));
        // d->array = xnp1->array;
        // memcpy(d->array, xnp1->array, xnp1->size * sizeof(unsigned long long));
        bignum *ac = mul(a, c);
        bignum *bd = mul(b, d);
        bignum *a_add_b = add(a, b);
        bignum *c_add_d = add(c, d);

        bignum *abcd = mul(a_add_b, c_add_d);
        free(a_add_b);
        free(c_add_d);
        bignum *abcd_sub_ac = sub(abcd, ac);
        bignum *ad_add_bc = sub(abcd_sub_ac, bd);
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
        bignum *res = add(add(ac, ad_add_bc), bd);
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

bignum *mul2num(elem_size_t x, elem_size_t y)
{
    if (x == 0 || y == 0)
    {
        bignum *res = malloc(sizeof(*res));
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
        free(ac_add_ad);
        bignum *res = add(ac_add_ad_add_bc, &bd);
        zero_justify(res);
        return res;
    }
}

/**
 *
 * Divide big number xn by big number xnp1.
 *
 */

bignum *div2bignums(bignum *xn, bignum *xnp1)
{
    size_t size = xn->size;
    bignum *res = new_bignum(size);
    bignum *row = new_bignum(size);
    row->size = 1;
    for (long i = xn->size - 1; i >= 0; i--)
    {
        array_shift(row, 1);
        row->array[0] = xn->array[i];
        res->array[i] = 0;
        while (compare_bignum(row, xnp1) < 1)
        {
            res->array[i]++;
            row = sub(row, xnp1);
            zero_justify(row);
        }
    }
    zero_justify(res);
    free(row);
    return res;
}

int compare_bignum(bignum *xn, bignum *xnp1)
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

void array_shift(bignum *n, int count)
{
    long i;
    elem_size_t *tmp = realloc(n, n->size + (count * sizeof(size_t)));
    ;
    if (tmp == NULL)
    {
        fprintf(stderr, "Couldn't reallocate memory for a n in arrayShift");
        exit(1);
    }
    else
    {
        n = tmp;
    }
    if ((n->size == 1) && (n->array[0] == 0))
        return;

    for (i = n->size - 1; i >= 0; i--)
        n->array[i + count] = n->array[i];

    for (i = 0; i < count; i++)
        *(n->array + i) = 0;
    n->size = n->size + count;
}

void zero_justify(bignum *n)
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
    const double log2_10 = 1.65;
    return (uint64_t)ceil(1 + log2_10 * numDigits);
}

/**
 *
 * Converts bignum to hexadecimal format string
 *
 */

char *hexToPrint(bignum *a)
{
    int count = a->size;
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
    if (!(str = string = malloc(len)))
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

/* void sqrt2(uint64_t n, bignum *xn, bignum *xnp1)
{
}*/
