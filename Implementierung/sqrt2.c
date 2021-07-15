#include "assert.h"
#include "errno.h"
#include "math.h"
#include "time.h"
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

extern void sum(uint64_t n, bignum *xn, bignum *xnp1);
bignum *mul(bignum *xn, bignum *xnp1);
void zeroJustify(bignum *n);
void arrayShift(bignum *n, int count);
int compareBignum(bignum *xn, bignum *xnp1);
bignum *div2bignums(bignum *xn, bignum *xnp1);
bignum *add(bignum *xn, bignum *xnp1);
bignum *sub(bignum *xn, bignum *xnp1);
bignum *mul2num(elem_size_t x, elem_size_t y);
matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2);
matrix *matrixBinaryExponentiation(unsigned long long n, int highestBit);
int calculateHighestBit(unsigned long long n);
uint64_t convertAccToN(uint64_t numDigits);
char *hexToPrint(bignum *a);
char *decToPrint(bignum *a);
void printUsage(char **argv);
void freeBigNum(bignum *toFree);
bignum *divideLongDivision(bignum *dividend, bignum *divisor);
bignum *copy(bignum *from, int countBits);
bignum *bitShiftRight(bignum *n, uint64_t count);
bignum *bitShiftLeft(bignum *n, uint64_t count);
bignum *addIntToBignum(bignum *n, int toAdd);
matrix *matrixSimpleExponentiation(unsigned long long n);
void freeMatrix(matrix *toFree);

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

int main(int argc, char **argv)
{
      int hexadecimal = 0;
      int c;
      if (argc == 1 || argc > 3)
      {
          printUsage(argv);
          return 1;
      }
      while (1)
      {
          int option_index = 0;
          c = getopt_long(argc, argv, "hx", long_options, &option_index);
          if (c == -1)
          {
              break;
          }
          switch (c)
          {
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
      uint64_t n = strtoull(argv[1], NULL, 0);
      if (errno == ERANGE)
      {
          fprintf(stderr, "the given number can not be represented, please pick a number < %llu", UINT64_MAX);
          return 1;
      }
      if (n == 0ULL || *argv[1] == '-')
      {
          fprintf(stderr, "invalid number of digits: it should be > 0 or the given parameter was not a number");
          return 1;
      }
      if (*argv[1] == '0') {
          fprintf(stderr, "Given number begins with 0, enter another number 0 < number < %llu", UINT64_MAX);
          return 1;
      }
    int op = convertAccToN(n);
    matrix *res = matrixBinaryExponentiation(op, calculateHighestBit(op));
    // double sqrt2 = 1 + (double)res->xn->array[0]/(double)res->xnp1->array[0];
    char *resultN = decToPrint(res->xn);
    char *resultNh = hexToPrint(res->xn);
    bignum *oper = new_bignum(1);
    *oper->array = 10;
    for (int i = 0; i < n; i++)
    {
        res->xn = mul(res->xn, oper);
    }
    // uint64_t s = res->xn->array[0] / res->xnp1->array[0];
    bignum *div = div2bignums(res->xn, res->xnp1);
    char *final = decToPrint(div);
    printf("%s\n", final);
    printf("size1: %u\n", res->xn->size);
    printf("xn: %s\n", resultN);
    printf("xnh: %s\n", resultNh);
    freeMatrix(res);
    // printf("sqrt2:%f\n", sqrt2);
    // printf("Result int res5->xn->array[0]: %u and res5-xn->array[1]: %u", res5->xn->array[0], res5->xn->array[1]);
    return 0;
}

/**
 *
 * Multiply two matrices.
 *
 */

void printUsage(char **argv)
{
    printf("Usage: %s <number of digits> \n", argv[0]);
    printf("-x: sets output numeral system to hexadecimal, default: decimal\n");
    printf("-h|--help: for usage help");
}

matrix *matrixMultiplication(matrix *matrix1, matrix *matrix2)
{
    matrix *res = malloc(3 * sizeof(struct bignum));
    if (res == 0)
    {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixMul");
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

matrix *matrixSimpleExponentiation(unsigned long long n)
{
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
    matrix->xnm1 = new_bignum(0);
    // matrix->xnm1->array[0] = 0;
    struct matrix *matrixInitial = malloc(3 * sizeof(bignum));
    if (matrixInitial == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for a matrix in matrixBinaryExp");
        exit(1);
    }
    matrixInitial->xnm1 = new_bignum(0);
    // *(matrixInitial->xnm1->array) = 0;
    matrixInitial->xn = new_bignum(1);
    *(matrixInitial->xn->array) = 1;
    matrixInitial->xnp1 = new_bignum(1);
    *(matrixInitial->xnp1->array) = 2;
    for (int i = 0; i < n - 1; i++)
    {
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

matrix *matrixBinaryExponentiation(unsigned long long n, int highestBit)
{
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
        struct matrix *tmp = matrixMultiplication(matrix, matrix);
        freeMatrix(matrix);
        matrix = tmp;
        if (n >> i & 1)
        {
            tmp = matrixMultiplication(matrix, matrixInitial);
            freeMatrix(matrix);
            matrix = tmp;
        }
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
    bignum *res = new_bignum(size2);
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
    zeroJustify(res);
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

    if (xn->size == 0 || xnp1->size == 0 || xn->size == 1 && *xn->array == 0 || xnp1->size == 1 && *xnp1->array == 0)
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
        uint64_t aNewSize;
        uint64_t bNewSize;
        uint64_t cNewSize;
        uint64_t dNewSize;
        if (size1 < size2)
        {
            aNewSize = size1 / 2;
            bNewSize = size1 - aNewSize;
            dNewSize = bNewSize;
            cNewSize = size2 - dNewSize;
        }
        else
        {
            cNewSize = size2 / 2;
            dNewSize = size2 - cNewSize;
            bNewSize = dNewSize;
            aNewSize = size1 - bNewSize;
        }
        bignum *b = new_bignum(bNewSize);

        // we need to write in b the second part of xn
        // void *memcpy(void *dest, const void *src, std::size_t count);
        // memcpy(b->array, xn->array, bNewSize * sizeof(elem_size_t));
        b->array = xn->array;

        bignum *a = new_bignum(aNewSize);
        // memcpy(a->array, (xn->array + bNewSize), aNewSize * sizeof(elem_size_t));
        a->array = xn->array + bNewSize;
        bignum *c = new_bignum(cNewSize);
        c->array = xnp1->array + dNewSize;
        // memcpy(c->array, (xnp1->array + dNewSize), cNewSize * sizeof(elem_size_t));

        bignum *d = new_bignum(dNewSize);
        // memcpy(d->array, xnp1->array, dNewSize * sizeof(elem_size_t));
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
        if (ac->size != 0)
        {
            arrayShift(ac, bNewSize + dNewSize);
            /* elem_size_t *tmp = realloc(ac->array, sizeof(elem_size_t) * (ac->size + 2));
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
            }*/
        }
        if (ad_add_bc->size != 0)
        {
            /* elem_size_t *tmp = realloc(ad_add_bc->array, sizeof(elem_size_t) * (ad_add_bc->size + 1));
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
            }*/
            arrayShift(ad_add_bc, bNewSize);
        }
        bignum *res = add(add(ac, ad_add_bc), bd);
        freeBigNum(ac);
        freeBigNum(bd);
        freeBigNum(abcd);
        freeBigNum(abcd_sub_ac);
        freeBigNum(ad_add_bc);
        zeroJustify(res);

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

bignum *div2bignums(bignum *xn, bignum *xnp1)
{
    bignum *res = new_bignum(0);
    res->array[0] = 0;
    uint64_t xnBits = xn->size >= 2 ? (xn->size - 1) * N + calculateHighestBit(xn->array[xn->size - 1]) + 1
                                    : calculateHighestBit(xn->array[0]) + 1;
    uint64_t xnp1Bits = xnp1->size >= 2 ? (xnp1->size - 1) * N + calculateHighestBit(xnp1->array[xnp1->size - 1] + 1)
                                        : calculateHighestBit(xnp1->array[0]) + 1;
    uint64_t diff = xnBits - xnp1Bits;
    bignum *temp2;
    while (diff >= 0) {
        bignum *toAdd = new_bignum(1);
        toAdd->array[0] = 1;
        bitShiftLeft(toAdd, diff);
        temp2 = mul(xnp1, toAdd);
        int compare = compareBignum(xn, temp2);
        if (compare < 1) {
            res = add(res, toAdd);
            xn = sub(xn, temp2);
            //freeBigNum(fitInitial);
            //fitInitial = new_bignum(1);
            if (compare == 0) {
                break;
            }
        }
        if (diff == 0) {
            break;
        }
        diff--;
    }
    free(temp2);
    /*while (diff > 1)
    {
        bignum *tmp = new_bignum(1);
        tmp->array[0] = 1;

        bitShiftLeft(tmp, diff - 1);
        res = add(res, tmp);
        bignum *july = mul(xnp1, tmp);
        xn = sub(xn, july);
        // july = mul(july, temp2);
        xnBits = xn->size >= 2 ? (xn->size - 2) * N + calculateHighestBit(xn->array[xn->size - 1])
                               : calculateHighestBit(xn->array[0]);
        xnp1Bits = xnp1->size >= 2 ? (xnp1->size - 2) * N + calculateHighestBit(xnp1->array[xnp1->size - 1])
                                   : calculateHighestBit(xnp1->array[0]);
        diff = xnBits - xnp1Bits;
    }
     */
    return res;
}

int compareBignum(bignum *xn, bignum *xnp1)
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

void arrayShift(bignum *n, int count)
{
    long i;
    elem_size_t *tmp = realloc(n->array, (n->size + count) * sizeof(elem_size_t));
    if (tmp == NULL)
    {
        fprintf(stderr, "Couldn't reallocate memory for a n in arrayShift");
        exit(1);
    }
    else
    {
        n->array = tmp;
    }
    if ((n->size == 0) && (n->array[0] == 0))
        return;

    for (i = n->size - 1; i >= 0; i--)
        n->array[i + count] = n->array[i];

    for (i = 0; i < count; i++)
        *(n->array + i) = 0;
    n->size = n->size + count;
    zeroJustify(n);
}

void zeroJustify(bignum *n)
{
    while ((n->size == 1 && n->array[0] == 0) || (n->size > 1) && (n->array[n->size - 1] == 0))
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
    return 2 + (LOG2_10 / 2) * numDigits;
}

/**
 *
 * Converts bignum to hexadecimal format string
 *
 */

char *hexToPrint(bignum *a)
{
    uint64_t count = a->size;
    if (count == 0)
    {
        char *zero = malloc(2);
        if (zero == NULL)
        {
            fprintf(stderr, "Couldn't allocate memory for a zerostring");
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
    if (tmp == 0)
    {
        fprintf(stderr, "Couldn't allocate memory for a tmp");
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
    if (string == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for a string in hexPrint");
        exit(1);
    }
    snprintf(str, addLen + 1, "%x", first);
    str += addLen;
    if (count == 0)
    {
        *str = '\0';
        return string;
    }
    count--;
    while (1)
    {
        elem_size_t c = *(array + count);
        snprintf(str, numCharinOneElem + 1, "%08x", c);
        str += numCharinOneElem;
        if (count == 0)
        {
            break;
        }
        else
        {
            count--;
        }
    }
    *str = '\0';
    return string;
}

void addToPrint(elem_size_t *decArray, uint64_t *size, elem_size_t summand)
{
    uint64_t tmp = summand;
    for (int i = 0; i < *size; i++)
    {
        if (tmp == 0)
        {
            return;
        }
        *(decArray + i) = (tmp += *(decArray + i)) % POWER10_9;
        tmp /= POWER10_9;
    }
    if (tmp)
    {
        *(decArray + *size) = tmp;
        (*size)++;
    }
}

void mulToPrint(elem_size_t *decArray, uint64_t *size)
{
    uint64_t tmp = 0;
    uint64_t max = ELEM_SIZE_MAX;
    max++;
    for (int i = 0; i < *size; i++)
    {
        *(decArray + i) = (tmp += *(decArray + i) * max) % POWER10_9;
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

char *decToPrint(bignum *a)
{
    uint64_t count = a->size;
    if (count == 0)
    {
        char *zero = malloc(2);
        if (zero == NULL)
        {
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
    if (decArray == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for decArray in decToPrint");
        exit(1);
    }
    elem_size_t *array = a->array;
    *decArray = 0;
    size = 1;
    uint64_t tmp = *(array + count);
    if (tmp != 0)
    {
        addToPrint(decArray, &size, tmp);
    }
    if (count > 0)
    {

        count--;
        while (1)
        {
            mulToPrint(decArray, &size);
            tmp = *(array + count);
            if (tmp)
            {
                addToPrint(decArray, &size, tmp);
            }
            if (count == 0)
            {
                break;
            }
            else
            {
                count--;
            }
        }
    }
    tmp = *(decArray + size - 1);
    elem_size_t temp = tmp;
    char *string, *str;
    char *toCheck = malloc(11);
    if (toCheck == 0)
    {
        fprintf(stderr, "Couldn't allocate memory for a tmp");
        exit(1);
    }
    // always add plus 1 byte because snprintf always adds and of the string, but we overwrite it then
    uint16_t addLen = snprintf(toCheck, 11, "%u", temp);
    free(toCheck);
    size--;
    uint64_t notOF = size;
    uint64_t max = 10 * notOF + addLen + 1;
    string = str = malloc(max);

    if (string == NULL)
    {
        fprintf(stderr, "Couldn't allocate memory for string in decToPrint");
        exit(1);
    }
    snprintf(str, addLen + 1, "%u", temp);
    str += addLen;
    if (size == 0)
    {
        *str = '\0';
        free(decArray);
        return string;
    }
    size--;
    while (1)
    {
        temp = *(decArray + size);
        str += snprintf(str, 11, "%09u", temp);
        if (size == 0)
        {
            break;
        }
        else
        {
            size--;
        }
    }
    *str = '\0';
    // free(decArray);
    return string;
}

/**
 *
 * Frees the memory
 *
 */

void freeBigNum(bignum *toFree)
{
    if (toFree != NULL)
    {
        if (toFree->size > 0)
        {
            free(toFree->array);
        }
        free(toFree);
    }
}

void freeMatrix(matrix *toFree)
{
    if (toFree != NULL)
    {
        freeBigNum(toFree->xn);
        freeBigNum(toFree->xnp1);
        freeBigNum(toFree->xnm1);
        free(toFree);
    }
}

bignum *divideLongDivision(bignum *dividend, bignum *divisor)
{
    int highestBitDividend = calculateHighestBit(dividend->array[dividend->size - 1]);
    int highestBitDivisor = calculateHighestBit(divisor->array[divisor->size - 1]);
    uint64_t k = N * (dividend->size - 1) + highestBitDividend + 1;
    uint64_t l = N * (divisor->size - 1) + highestBitDivisor + 1;
    bignum *d = new_bignum(divisor->size);
    bignum *r = copy(dividend, l - 1);
    bignum *q = new_bignum(dividend->size);
    bignum *arrayTemp = new_bignum(r->size);
    for (int i = 0; i <= k - l; i++)
    {
        uint64_t arrayNumber = dividend->size - 1 - ((31 - highestBitDividend + i + l - 1) / N);
        int64_t bitInArray;
        if (arrayNumber == dividend->size - 1)
        {

            bitInArray = highestBitDividend - (i + l - 1) % N;
        }
        else
        {
            bitInArray = 32 - (i + l - 1) % N;
        }
        /* if (bitInArray < 0)
        {
            bit
            InArray += N;
        }
        */
        int alpha = dividend->array[arrayNumber] & (1 << bitInArray);
        memcpy(arrayTemp->array, r->array, (r->size + 1) * sizeof(elem_size_t));
        arrayTemp = bitShiftLeft(r, 1);
        d->array = arrayTemp->array;
        if (alpha > 0)
        {
            d = addIntToBignum(d, 1);
        }
        int subtract = compareBignum(d, divisor);
        if (subtract < 1)
        {
            r = sub(d, divisor);
        }
        else
        {
            memcpy(r->array, d->array, (d->size) * sizeof(elem_size_t));
        }
        bitShiftLeft(q, 1);
        if (subtract < 1)
        {
            q = addIntToBignum(q, 1);
        }
    }
    printf("Quotient: %u and rest: %u\n", q->array[0], r->array[0]);
    free(arrayTemp);
    free(d);
    return q;
}

bignum *bitShiftRight(bignum *n, uint64_t count)
{
    if (n->size == 0 && n->array[0] == 0)
    {
        return n;
    }
    while (count > 0)
    {
        for (int i = 0; i < n->size; i++)
        {
            n->array[i] = n->array[i] >> 1;
            if (i != n->size - 1 && n->array[i + 1] & 1)
            {
                n->array[i] |= 1 << (N - 1);
            }
        }
        count--;
    }
    zeroJustify(n);
    return n;
}

bignum *bitShiftLeft(bignum *n, uint64_t count)
{
    if (n->size == 0 && n->array[0] == 0)
    {
        return n;
    }
    while (count > 0)
    {
        if (count > 31)
        {
            count -= 32;
            arrayShift(n, 1);
        }
        else
        {
            for (int i = n->size - 1; i >= 0; i--)
            {
                if (n->array[i] >= TWO_POW_31)
                {
                    if (i == n->size - 1)
                    {
                        n->array = realloc(n->array, (n->size + 1) * sizeof(elem_size_t));
                        n->size++;
                        n->array[n->size - 1] = 0;
                    }
                    n->array[i + 1] += 1;
                }
                n->array[i] = n->array[i] << 1;
            }
            count--;
        }
    }
    return n;
}

bignum *addIntToBignum(bignum *n, int toAdd)
{
    bignum *temp = new_bignum(1);
    temp->array[0] = toAdd;
    return add(n, temp);
}

bignum *copy(bignum *from, int countBits)
{
    int size = (countBits - 1) / N;
    bignum *res = new_bignum(size + 1);
    int currentArrayFrom = from->size - 1;
    int currentBitFrom = calculateHighestBit(from->array[from->size - 1]);
    int currentArraySet = (countBits - 1) / N;
    int currentBitSet = (countBits - 1) % N;
    for (int i = countBits - 1; i >= 0; i--)
    {
        if (from->array[currentArrayFrom] & 1ULL << currentBitFrom)
        {
            res->array[currentArraySet] |= 1ULL << currentBitSet;
        }
        if (currentBitFrom == 0)
        {
            currentArrayFrom--;
            currentBitFrom = 31;
        }
        else
        {
            currentBitFrom--;
        }
        if (currentBitSet == 0)
        {
            currentArraySet--;
            currentBitSet = 31;
        }
        else
        {
            currentBitSet--;
        }
    }
    return res;
}