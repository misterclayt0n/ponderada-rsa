/*
 * Test cases for Trial Division and Pollard's Rho factorization algorithms
 * Usage: ./test_factorization
 */

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

// ============ Trial Division ============
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

// ============ Pollard's Rho ============
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

uint64_t f(uint64_t x, uint64_t n)
{
    return ((__uint128_t)x * x + 1) % n;
}

uint64_t pollards_rho(uint64_t n, uint64_t *iterations)
{
    *iterations = 0;
    
    if (n % 2 == 0)
    {
        *iterations = 1;
        return 2;
    }
    
    uint64_t x = 2, y = 2, d = 1;
    
    while (d == 1)
    {
        (*iterations)++;
        x = f(x, n);
        y = f(f(y, n), n);
        
        uint64_t diff = (x > y) ? x - y : y - x;
        d = gcd(diff, n);
        
        if (*iterations > 10000000)
            return 0;
    }
    
    return (d != n) ? d : 0;
}

// ============ Test Framework ============
typedef struct {
    uint64_t n;
    uint64_t expected_p;
    uint64_t expected_q;
    const char *description;
} TestCase;

int test_algorithm(const char *name, uint64_t (*factor_func)(uint64_t, uint64_t*), TestCase *tests, int num_tests)
{
    int passed = 0;
    int failed = 0;
    
    printf("Testing %s\n", name);
    printf("----------------------------------------\n");
    
    for (int i = 0; i < num_tests; i++)
    {
        uint64_t iterations;
        uint64_t p = factor_func(tests[i].n, &iterations);
        uint64_t q = (p != 0 && p != tests[i].n) ? tests[i].n / p : 0;
        
        // Check if we found valid factors (order doesn't matter)
        int correct = 0;
        if (p != 0 && q != 0)
        {
            if ((p == tests[i].expected_p && q == tests[i].expected_q) ||
                (p == tests[i].expected_q && q == tests[i].expected_p))
            {
                correct = 1;
            }
            // Also accept if p * q == n (valid factorization even if different factors)
            if (p * q == tests[i].n && p > 1 && q > 1)
            {
                correct = 1;
            }
        }
        
        if (correct)
        {
            printf("  [PASS] %s: %" PRIu64 " = %" PRIu64 " * %" PRIu64 "\n", 
                   tests[i].description, tests[i].n, p, q);
            passed++;
        }
        else
        {
            printf("  [FAIL] %s: %" PRIu64 "\n", tests[i].description, tests[i].n);
            printf("         Expected: %" PRIu64 " * %" PRIu64 "\n", 
                   tests[i].expected_p, tests[i].expected_q);
            printf("         Got: %" PRIu64 " * %" PRIu64 "\n", p, q);
            failed++;
        }
    }
    
    printf("----------------------------------------\n");
    printf("Results: %d passed, %d failed\n\n", passed, failed);
    
    return failed;
}

int main()
{
    printf("Factorization Algorithm Test Suite\n");
    printf("========================================\n\n");
    
    TestCase tests[] = {
        // Small semiprimes
        {15, 3, 5, "Small: 3 * 5"},
        {35, 5, 7, "Small: 5 * 7"},
        {77, 7, 11, "Small: 7 * 11"},
        {91, 7, 13, "Small: 7 * 13"},
        {143, 11, 13, "Small: 11 * 13"},
        {221, 13, 17, "Small: 13 * 17"},
        
        // Medium semiprimes
        {3233, 53, 61, "Medium: 53 * 61"},
        {5767, 73, 79, "Medium: 73 * 79"},
        {10403, 101, 103, "Medium: 101 * 103"},
        {19043, 137, 139, "Medium: 137 * 139"},
        
        // Larger semiprimes
        {129834181, 5573, 23297, "Large: 5573 * 23297"},
        {1106774983, 32771, 33773, "Large: 32771 * 33773"},
        {3215031751ULL, 56711, 56701, "Large: 56711 * 56701"},
        
        // Even larger (64-bit safe)
        {275447306077ULL, 524309, 525353, "XLarge: 524309 * 525353"},
        {4400626126189ULL, 2097257, 2098277, "XLarge: 2097257 * 2098277"},
        
        // Edge cases
        {4, 2, 2, "Edge: 2 * 2"},
        {6, 2, 3, "Edge: 2 * 3"},
        {9, 3, 3, "Edge: 3 * 3"},
        {49, 7, 7, "Edge: 7 * 7"},
        
        // RSA-like (balanced primes)
        {70377803883943ULL, 8388617, 8389679, "RSA-like: balanced 24-bit primes"},
    };
    
    int num_tests = sizeof(tests) / sizeof(tests[0]);
    
    int td_failures = test_algorithm("Trial Division", trial_division, tests, num_tests);
    int pr_failures = test_algorithm("Pollard's Rho", pollards_rho, tests, num_tests);
    
    printf("========================================\n");
    printf("Final Summary\n");
    printf("========================================\n");
    printf("Trial Division: %d/%d tests passed\n", num_tests - td_failures, num_tests);
    printf("Pollard's Rho:  %d/%d tests passed\n", num_tests - pr_failures, num_tests);
    printf("\n");
    
    if (td_failures == 0 && pr_failures == 0)
    {
        printf("All tests passed!\n");
        return 0;
    }
    else
    {
        printf("Some tests failed.\n");
        return 1;
    }
}
