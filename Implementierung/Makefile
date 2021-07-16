.PHONY: all clean

CC:=gcc
CFLAGS+=-Wall -Wextra -std=c11 -Wpedantic -D_POSIX_SOURCE -D_GNU_SOURCE

sqrt2: sqrt2.c
	gcc ${CFLAGS} -o $@ $^ -lm

all: sqrt2

clean: rm -rf main
