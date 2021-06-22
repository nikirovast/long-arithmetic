.intel_syntax noprefix
.global quicksort_asm
.text

// void sqrt2(uint64_t n, struct bignum *xn, struct bignum *xnp1);
// This function adds two numbers assumming, there's enough space and the first one is bigger or has the same size as the second
// The result is in the first parameter
// !IMP First param should have enough space
// Could be optimized later
sum:
    mov r12, [r8]
    mov r13, [rdx]
    mov rax, 0 
    summ:
    cmp rax, r12
    mov rbx, 8
    je loop_end
    add rdx, rbx
    add r8, rbx
    mov r14, [rdx]
    mov r15, [r8]
    add r14, rsi
    add r14, r15
    mov rsi, 0
    adc rsi, 0
    mov[rdx], r14
    inc rax
    jmp summ
    loop_end:
    cmp r12, r13
    je sum_end
    add rdx, rbx
    mov r14, [rdx]
    add r14, rsi
    mov[rdx], r14
    sum_end:
    ret