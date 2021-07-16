.PHONY: all clean

CC:=gcc
CFLAGS+=-Wall -Wextra -std=c11 -Wpedantic -D_POSIX_SOURCE -D_GNU_SOURCE

main: sqrt2.c
	gcc ${CFLAGS} -o $@ $^ -lm

all: main

clean: rm -rf main