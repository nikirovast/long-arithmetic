//
// Created by Andrey Manucharyan on 20.06.21.
//
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errno.h"
#include "assert.h"

#define N 32
#define ELEM_SIZE_MAX UINT32_MAX
typedef uint32_t elem_size_t;
typedef struct bignum bignum;

struct bignum {
  size_t size;
  elem_size_t* array;
};

extern void sum(uint64_t n, bignum* xn, bignum* xnp1);
bignum* mul(bignum* xn, bignum* xnp1);
// void sqrt2(uint64_t n, bignum *xn, bignum *xnp1);#
bignum* add(bignum* xn, bignum* xnp1);
bignum* sub(bignum* xn, bignum* xnp1);
bignum* mul2num(elem_size_t x, elem_size_t y);

bignum* new_bignum(size_t size) {
  bignum* res = malloc(sizeof(bignum));
  if (res == NULL) {
    fprintf(stderr, "Couldn't allocate memory for a struct");
    exit(1);
  }
  res->size = size;
  res->array = malloc(sizeof(elem_size_t) * size);
  if (res->array == NULL) {
    fprintf(stderr, "Couldn't allocate memory for an array");
    exit(1);
  }
  return res;
}

// in VSCode 1st is in RCX, 2nd in RDX, 3rd in R8
int main(int argc, char** argv) {
  // Yulia for tests

  bignum *a = new_bignum(2);
  bignum *b = new_bignum(1);
  *(a->array + 0) = 4294967295;
  *(a->array + 1) = 1;
  *(b->array + 0) = 4294967295;
    printf("%u, %u, %u", a->array[0], a->array[1], b-> array[0]);
  bignum *res = mul(a, b);

  uint64_t s = res->size;

    printf("\n");
    printf("Size: %lu\n",res->size);
  for (int i = res-> size - 1; i >= 0; i--) {
      printf("%u\n", res->array[i]);
  }


  //test 1
    bignum *first = new_bignum(2);
    bignum *second = new_bignum(1);
    *(first->array + 0) = 1234;
    *(first->array + 1) = 1;
    *(second->array + 0) = 129;
    bignum *res1 = mul(first, second);
    //8589937060
    assert(res1->array[0] == 2468);
    assert(res1->array[1] == 2);

    //test 2
    bignum *third = new_bignum(1);
    bignum *fourth = new_bignum(1);
    *(third->array) = 65537;
    *(fourth->array) = 89340;
    bignum *res2 = mul(third, fourth);
    //5855075580
    assert(res2->array[0] == 1560108285);
    assert(res2->array[1] == 1);

  /* bignum *b = new_bignum(sizeof(elem_size_t) * numItems);
  b->size = numItems;
  uint64_t arr2[2] = {3, 3};
  strcpy(a->array, arr2);
  sum(3, a, b);
   */

  // implement parsing of arguments, size = 5 as example
  if (argc != 3) {
    fprintf(stderr,
            "usage: %s <number of iterations> <decimal (d) or hexadecimal (x) "
            "output>\n",
            argv[0]);
    return 1;
  }
  uint64_t size = strtoull(argv[1], NULL, 0);
  if (errno == ERANGE) {
    fprintf(stderr,
            "the given number can not be represented, please pick a number < "
            "UINT64_MAX");
    return 1;
  }
  if (size == 0ULL || *argv[1] == '-') {
    fprintf(stderr,
            "invalid number of iterations: it should be > 0 or the given "
            "parameter was not a number");
    return 1;
  }
  char* output = argv[2];
  if ((*output != 'd' && *output != 'x') || strlen(output) > 1) {
    fprintf(stderr, "invalid numeral system for output");
    return 1;
  }
  bignum* xn =
      new_bignum(size + 1);  // size + 1 because of size parameter in the struct
  if (xn == NULL) {
    fprintf(stderr,
            "not enough memory for the given number of iterations: %llu", size);
    return 1;
  }
  xn->array[0] = 1;

  bignum* xnp1 =
      new_bignum(size + 1);  // size + 1 because of size parameter in the struct
  if (xnp1 == NULL) {
    fprintf(stderr,
            "not enough memory for the given number of iterations: %llu", size);
    return 1;
  }
  xnp1->array[0] = 2;
  printf(
      "memory for xn and xnp1 successfully allocated, value in the first "
      "element of array with size %llu: xn is "
      "%llu and xnp1 is %llu.\n",
      size, xn->array[0], xnp1->array[0]);
  if (*output == 'd') {
    printf("output numeral system is decimal");
  } else {
    printf("output numeral system is hexadecimal");
  }
  return 0;
}

/**
 *
 * Sum two big numbers.
 *
 */

bignum* add(bignum* xn, bignum* xnp1) {
  if (xn->size > xnp1->size) {
    bignum* tmp = xn;
    xn = xnp1;
    xnp1 = tmp;
  }

  uint64_t size1 = xn->size;
  uint64_t size2 = xnp1->size;
  bignum* res = new_bignum(sizeof(elem_size_t) * size2);
  res->size = size2;
  int i;
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
    elem_size_t* tmp = realloc(res->array, sizeof(elem_size_t) * (size2 + 1));
    if (tmp == NULL) {
      fprintf(stderr, "Couldn't reallocate memory for a carry in sum");
      exit(1);

    } else {
      res->array = tmp;
      free(tmp);
      res->size += 1;
      *(res->array + i) = carry;
    }
  }
  return res;
}

/**
 *
 * Subtract xnp1 from xn.
 *
 */

bignum* sub(bignum* xn, bignum* xnp1) {
  uint64_t size1 = xn->size;
  uint64_t size2 = xnp1->size;
  bignum* res = new_bignum(sizeof(elem_size_t) * size1);
  res->size = size1;
  int i;
  elem_size_t carry = 0;
  for (i = 0; i < size2; i++) {
    elem_size_t a = *(xn->array + i);
    elem_size_t b = *(xnp1->array + i);
    if (a < 1 && (b > 0 || carry > 0) || a - carry < b) {
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
  return res;
}

/**
 *
 * Multiply two big numbers xnp1 by xn.
 *
 */

bignum* mul(bignum* xn, bignum* xnp1) {
  uint64_t size1 = xn->size;
  uint64_t size2 = xnp1->size;

  if (size1 > 1 || size2 > 1) {
    uint64_t aNewSize = size1 / 2;
    bignum* a = malloc(sizeof(*a));
    if (a == NULL) {
      fprintf(stderr, "Couldn't allocate memory for a in mul");
      exit(1);
    }
    a->size = aNewSize;
    a->array = xn->array;
    uint64_t bNewSize = size1 - aNewSize;
    bignum* b = malloc(sizeof(*b));
    if (b == NULL) {
      fprintf(stderr, "Couldn't allocate memory for b in mul");
      exit(1);
    }
    b->size = bNewSize;
    b->array = xn->array + aNewSize;

    uint64_t cNewSize = size2 / 2;
    bignum* c = malloc(sizeof(*c));
    if (c == NULL) {
      fprintf(stderr, "Couldn't allocate memory for c in mul");
      exit(1);
    }
    c->size = cNewSize;
    c->array = xnp1->array;
    uint64_t dNewSize = size2 - cNewSize;
    bignum* d = malloc(sizeof(*d));
    if (d == NULL) {
      fprintf(stderr, "Couldn't allocate memory for d in mul");
      exit(1);
    }
    d->size = dNewSize;
    d->array = xnp1->array + dNewSize;

    bignum* ac = mul(a, c);
    bignum* bd = mul(b, d);
    bignum* a_add_b = add(a, b);
    bignum* c_add_d = add(c, d);
    free(a);
    free(b);
    free(c);
    free(d);
    bignum* abcd = mul(a_add_b, c_add_d);
    free(a_add_b);
    free(c_add_d);
    bignum* abcd_sub_ac = sub(abcd, ac);
    bignum* ad_add_bc = sub(abcd_sub_ac, bd);

    // TODO rewrite it here
     /*bignum* ac_shift =
    realloc(ac, sizeof(*ac) + sizeof(elem_size_t) * (ac->size + 2));*/

    elem_size_t* tmp = realloc(ac->array, sizeof(elem_size_t) * (ac->size + 2));
    if (tmp == NULL) {
      fprintf(stderr, "Couldn't reallocate memory for a shift in sum");
      exit(1);
    } else {
      ac->array = tmp;
      free(tmp);
      ac->size += 2;
    }
    tmp =
        realloc(ad_add_bc->array, sizeof(elem_size_t) * (ad_add_bc->size + 1));
    if (tmp == NULL) {
      fprintf(stderr, "Couldn't reallocate memory for a shift in sum");
      exit(1);
    } else {
      ad_add_bc->array = tmp;
      free(tmp);
      ad_add_bc->size += 1;
    }
     /*bignum* ad_add_bc_shift =
        realloc(ad_add_bc,
            sizeof(*ad_add_bc) + sizeof(elem_size_t) * (ad_add_bc->size + 1));*/
    bignum* res = add(add(ac_shift, ad_add_bc_shift), bd);
    free(abcd);
    free(abcd_sub_ac);
    free(ad_add_bc);
    return res;
  } else if (xn->size == 0 || xnp1->size == 0) {
    bignum* res = malloc(sizeof(*xn));
    res->size = 0;
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

bignum* mul2num(elem_size_t x, elem_size_t y) {
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
  elem_size_t ad_array[] = {(a * d << shift_size) >> shift_size,
                            (a * d) >> shift_size};
  ad.array = ad_array;

  bignum bc;
  bc.size = 2;
  elem_size_t bc_array[] = {(b * c << shift_size) >> shift_size,
                            (b * c) >> shift_size};
  bc.array = bc_array;

  bignum bd;
  bd.size = 1;
  elem_size_t bd_array[] = {b * d};
  bd.array = bd_array;
  bignum* ac_add_ad = add(&ac, &ad);
  bignum* ac_add_ad_add_bc = add(ac_add_ad, &bc);
  free(ac_add_ad);
  bignum* res = add(ac_add_ad_add_bc, &bd);
  free(ac_add_ad_add_bc);

  return res;
}
/* void sqrt2(uint64_t n, bignum *xn, bignum *xnp1)
{
}*/
