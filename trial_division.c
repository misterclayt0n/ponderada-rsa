/*
 * Trial Division Attack on RSA
 * Usage: ./trial_division <n> [e]
 *        ./trial_division --demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <time.h>

uint64_t gcd(uint64_t a, uint64_t b)
{
    while (b != 0)
    {
        uint64_t temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

int64_t mod_inverse(int64_t e, int64_t phi)
{
    int64_t t = 0, newt = 1;
    int64_t r = phi, newr = e;
    
    while (newr != 0)
    {
        int64_t quotient = r / newr;
        int64_t temp_t = t;
        t = newt;
        newt = temp_t - quotient * newt;
        int64_t temp_r = r;
        r = newr;
        newr = temp_r - quotient * newr;
    }
    
    if (t < 0)
        t += phi;
    return t;
}

uint64_t trial_division(uint64_t n, uint64_t *iterations)
{
    *iterations = 0;
    
    if (n % 2 == 0)
    {
        *iterations = 1;
        return 2;
    }
    
    uint64_t limit = (uint64_t)sqrt((double)n) + 1;
    for (uint64_t i = 3; i <= limit; i += 2)
    {
        (*iterations)++;
        if (n % i == 0)
            return i;
    }
    return n;
}

int is_prime(uint64_t n)
{
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    for (uint64_t i = 3; i * i <= n; i += 2)
        if (n % i == 0) return 0;
    return 1;
}

uint64_t next_prime(uint64_t n)
{
    while (!is_prime(n))
        n++;
    return n;
}

void run_demo()
{
    printf("Trial Division Scaling Demo\n");
    printf("============================\n\n");
    printf("%-10s %15s %12s %15s\n", "Bits", "Iterations", "Time", "Est. 1024-bit");
    printf("--------------------------------------------------------------\n");
    
    // Pre-computed n values with balanced primes valid for e=3
    struct { int bits; uint64_t n; } tests[] = {
        {16, 1106774983ULL},
        {20, 275447306077ULL},
        {22, 4400626126189ULL},
        {24, 70377803883943ULL},
        {26, 1125938964277027ULL},
        {28, 18014546685901351ULL},
        {30, 288230981742142951ULL},
        {31, 1152922614855900181ULL},
    };
    int num_tests = sizeof(tests) / sizeof(tests[0]);
    
    double last_time = 0;
    
    for (int i = 0; i < num_tests; i++)
    {
        clock_t start = clock();
        uint64_t iterations;
        uint64_t p = trial_division(tests[i].n, &iterations);
        clock_t end = clock();
        double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        
        // Estimate time for 1024-bit primes (2048-bit n)
        // Each bit doubles iterations, so 2^(1024 - bits) times slower
        int bits_remaining = 1024 - tests[i].bits;
        double est_seconds = time_spent * pow(2, bits_remaining);
        
        // Convert to years with all zeros
        double years = est_seconds / (365.25 * 24 * 3600);
        int exponent = (years > 0) ? (int)floor(log10(years)) : 0;
        
        printf("%-10d %15" PRIu64 " %10.4fs       ", tests[i].bits, iterations, time_spent);
        
        if (years < 1)
            printf("%.2f sec\n", est_seconds);
        else if (exponent < 10)
            printf("%.0f years\n", years);
        else
        {
            printf("1");
            for (int z = 0; z < exponent; z++)
                printf("0");
            printf(" years\n");
        }
        
        last_time = time_spent;
    }
    
    printf("\n");
    printf("Note: Real RSA uses 1024-bit primes (2048-bit n)\n");
    printf("Trial division is completely infeasible at that scale.\n");
    printf("\nUniverse age: ~13.8 billion years\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Usage: %s <n> [e]\n", argv[0]);
        printf("       %s --demo    (run scaling demonstration)\n", argv[0]);
        return 1;
    }
    
    if (strcmp(argv[1], "--demo") == 0)
    {
        run_demo();
        return 0;
    }
    
    uint64_t n = strtoull(argv[1], NULL, 10);
    uint64_t e = (argc >= 3) ? strtoull(argv[2], NULL, 10) : 3;
    
    if (n < 4)
    {
        fprintf(stderr, "Error: n must be >= 4\n");
        return 1;
    }
    
    printf("Trial Division Attack\n");
    printf("n = %" PRIu64 ", e = %" PRIu64 "\n\n", n, e);
    
    clock_t start = clock();
    uint64_t iterations;
    uint64_t p = trial_division(n, &iterations);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    
    if (p == n)
    {
        printf("Failed: n is prime\n");
        return 1;
    }
    
    uint64_t q = n / p;
    uint64_t phi = (p - 1) * (q - 1);
    
    printf("Factors: p = %" PRIu64 ", q = %" PRIu64 "\n", p, q);
    printf("Iterations: %" PRIu64 ", Time: %.6fs\n\n", iterations, time_spent);
    
    if (gcd(e, phi) != 1)
    {
        printf("Error: e is not valid for these primes\n");
        return 1;
    }
    
    int64_t d = mod_inverse(e, phi);
    
    printf("phi(n) = %" PRIu64 "\n", phi);
    printf("Private key d = %" PRId64 "\n\n", d);
    
    printf("Public:  (n=%" PRIu64 ", e=%" PRIu64 ")\n", n, e);
    printf("Private: (n=%" PRIu64 ", d=%" PRId64 ")\n", n, d);
    
    return 0;
}
