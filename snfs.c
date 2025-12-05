/*
 * Toy Special Number Field Sieve (SNFS) factorization
 * Usage:
 *   ./snfs <n> [e] [degree] [B] [K]
 *   ./snfs --demo
 *
 * Focus: educational, small semiprimes with special form n ~= m^degree + 1.
 * Defaults: degree=8, B=200 (factor base bound), K=5000 (search bound for k in 1-D sieve).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef __int128 i128;

// ============ Small helpers ============

u128 gcd_u128(u128 a, u128 b)
{
    while (b != 0)
    {
        u128 t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Prevent overflow by using double-and-add
u128 mul_mod(u128 a, u128 b, u128 mod)
{
    u128 res = 0;
    a %= mod;
    while (b)
    {
        if (b & 1)
        {
            res += a;
            if (res >= mod)
                res -= mod;
        }
        a <<= 1;
        if (a >= mod)
            a -= mod;
        b >>= 1;
    }
    return res % mod;
}

u128 pow_mod(u128 base, u128 exp, u128 mod)
{
    u128 result = 1;
    base %= mod;
    while (exp)
    {
        if (exp & 1)
            result = mul_mod(result, base, mod);
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

u128 pow_u128(u128 base, int exp)
{
    u128 res = 1;
    for (int i = 0; i < exp; i++)
        res *= base;
    return res;
}

// Parse decimal string to u128
u128 parse_u128(const char *s)
{
    u128 v = 0;
    for (int i = 0; s[i]; i++)
    {
        if (s[i] < '0' || s[i] > '9')
            continue;
        v = v * 10 + (s[i] - '0');
    }
    return v;
}

void print_u128(u128 x)
{
    char buf[64];
    int idx = 0;
    if (x == 0)
    {
        printf("0");
        return;
    }
    while (x > 0)
    {
        int digit = x % 10;
        buf[idx++] = '0' + digit;
        x /= 10;
    }
    for (int i = idx - 1; i >= 0; i--)
        putchar(buf[i]);
}

// Extended GCD for modular inverse (128-bit)
u128 mod_inverse_u128(u128 e, u128 phi)
{
    i128 t = 0, newt = 1;
    i128 r = (i128)phi, newr = (i128)e;
    
    while (newr != 0)
    {
        i128 q = r / newr;
        i128 temp_t = t;
        t = newt;
        newt = temp_t - q * newt;
        i128 temp_r = r;
        r = newr;
        newr = temp_r - q * newr;
    }
    if (t < 0)
        t += (i128)phi;
    return (u128)t;
}

// Integer d-th root: largest x with x^d <= n
u128 int_root(u128 n, int d)
{
    u128 low = 1, high = n, ans = 1;
    while (low <= high)
    {
        u128 mid = low + ((high - low) >> 1);
        u128 p = pow_u128(mid, d);
        if (p == 0 || p > n)
        {
            high = mid - 1;
        }
        else
        {
            ans = mid;
            low = mid + 1;
        }
    }
    return ans;
}

// ============ Prime generation ============

#define MAX_FB 6000   // max primes in factor base (primes <= ~60000)
#define LP_BOUND 100000000

int generate_primes(int limit, uint32_t primes[MAX_FB])
{
    int count = 0;
    char *is_prime = calloc(limit + 1, 1);
    if (!is_prime)
        return 0;
    for (int i = 2; i <= limit; i++)
        is_prime[i] = 1;
    for (int p = 2; p * p <= limit; p++)
    {
        if (is_prime[p])
        {
            for (int j = p * p; j <= limit; j += p)
                is_prime[j] = 0;
        }
    }
    for (int i = 2; i <= limit && count < MAX_FB; i++)
    {
        if (is_prime[i])
            primes[count++] = (uint32_t)i;
    }
    free(is_prime);
    return count;
}

// ============ Relation / matrix handling ============

#define MAX_REL 12000
#define MAX_EXP 8   // exponent counters stored in uint8_t

typedef struct {
    int a_offset;                // k such that a = m + k
    uint8_t r_exp[MAX_FB];       // exponents on rational side
    uint8_t a_exp[MAX_FB];       // exponents on algebraic side
} Relation;

static Relation relations[MAX_REL];
static int relation_count = 0;

// Bit matrix helpers
static uint64_t row_bits[MAX_REL][(2 * MAX_FB + 63) / 64];
static uint64_t combo_bits[MAX_REL][(MAX_REL + 63) / 64];
static int pivot_col[MAX_REL];
static int matrix_rows = 0;

static int first_set_bit(uint64_t *row, int words)
{
    for (int w = 0; w < words; w++)
    {
        if (row[w])
        {
            int offset = __builtin_ctzll(row[w]);
            return w * 64 + offset;
        }
    }
    return -1;
}

static void xor_rows(uint64_t *dst, uint64_t *src, int words)
{
    for (int i = 0; i < words; i++)
        dst[i] ^= src[i];
}

static int row_is_zero(uint64_t *row, int words)
{
    for (int i = 0; i < words; i++)
        if (row[i])
            return 0;
    return 1;
}

// Attempt to insert row; if dependent, returns 1 and fills dependency mask
static int insert_row(uint64_t *row, uint64_t *combo, int col_words, int combo_words, uint64_t *out_dep)
{
    for (int r = 0; r < matrix_rows; r++)
    {
        int pc = pivot_col[r];
        int word = pc / 64;
        int bit = pc % 64;
        if (row[word] & ((uint64_t)1 << bit))
        {
            xor_rows(row, row_bits[r], col_words);
            xor_rows(combo, combo_bits[r], combo_words);
        }
    }
    if (row_is_zero(row, col_words))
    {
        memcpy(out_dep, combo, combo_words * sizeof(uint64_t));
        return 1; // dependency found
    }
    int pc = first_set_bit(row, col_words);
    if (pc < 0)
        return 0;
    memcpy(row_bits[matrix_rows], row, col_words * sizeof(uint64_t));
    memcpy(combo_bits[matrix_rows], combo, combo_words * sizeof(uint64_t));
    pivot_col[matrix_rows] = pc;
    matrix_rows++;
    return 0;
}

// ============ SNFS core ============

// Factor a value using the factor base; fill exp counters; return 1 if fully smooth
static int is_prime_u64(uint64_t x)
{
    if (x < 2) return 0;
    if (x % 2 == 0) return x == 2;
    for (uint64_t i = 3; i * i <= x; i += 2)
    {
        if (x % i == 0)
            return 0;
    }
    return 1;
}

static int factor_with_fb(u128 value, uint32_t *primes, int *fb_size, uint8_t *exp_out)
{
    for (int i = 0; i < *fb_size; i++)
    {
        uint32_t p = primes[i];
        while ((value % p) == 0)
        {
            value /= p;
            if (exp_out[i] < 250)
                exp_out[i]++; // keep small
        }
    }
    if (value == 1)
        return 1;
    
    // Large-prime variant (single extra prime <= LP_BOUND)
    if (value <= LP_BOUND && *fb_size < MAX_FB && is_prime_u64((uint64_t)value))
    {
        primes[*fb_size] = (uint32_t)value;
        exp_out[*fb_size] = 1;
        (*fb_size)++;
        return 1;
    }
    return 0;
}

// Build dependency -> compute square congruence
static u128 attempt_dependency(uint64_t *dep_mask, int dep_words, uint32_t *primes, int fb_size, u128 n)
{
    uint32_t total_r[MAX_FB] = {0};
    uint32_t total_a[MAX_FB] = {0};
    
    for (int i = 0; i < relation_count; i++)
    {
        int word = i / 64;
        int bit = i % 64;
        if (!(dep_mask[word] & ((uint64_t)1 << bit)))
            continue;
        for (int j = 0; j < fb_size; j++)
        {
            total_r[j] += relations[i].r_exp[j];
            total_a[j] += relations[i].a_exp[j];
        }
    }
    
    u128 x = 1;
    u128 y = 1;
    for (int j = 0; j < fb_size; j++)
    {
        if (total_r[j])
        {
            u128 exp = total_r[j] / 2;
            x = mul_mod(x, pow_mod(primes[j], exp, n), n);
        }
        if (total_a[j])
        {
            u128 exp = total_a[j] / 2;
            y = mul_mod(y, pow_mod(primes[j], exp, n), n);
        }
    }
    
    u128 diff = (x > y) ? (x - y) : (y - x);
    u128 g = gcd_u128(diff, n);
    if (g > 1 && g < n)
        return g;
    
    u128 sum = x + y;
    if (sum >= n)
        sum -= n;
    g = gcd_u128(sum, n);
    if (g > 1 && g < n)
        return g;
    
    return 0;
}

static uint64_t gcd_u64(uint64_t a, uint64_t b)
{
    while (b != 0)
    {
        uint64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// ============ Fallback: simple Pollard rho for u128 (educational only) ============
static u128 rho_func(u128 x, u128 c, u128 n)
{
    return (mul_mod(x, x, n) + c) % n;
}

static u128 pollard_rho_u128(u128 n)
{
    if ((n & 1) == 0)
        return 2;
    u128 c = 1;
    for (int attempt = 0; attempt < 5; attempt++, c += 2)
    {
        u128 x = 2, y = 2, d = 1;
        for (int i = 0; i < 200000; i++)
        {
            x = rho_func(x, c, n);
            y = rho_func(rho_func(y, c, n), c, n);
            u128 diff = (x > y) ? (x - y) : (y - x);
            d = gcd_u128(diff, n);
            if (d > 1 && d < n)
                return d;
        }
    }
    return 0;
}

u128 snfs_factor(u128 n, int degree, int fb_bound, int window)
{
    uint32_t primes[MAX_FB];
    int fb_size = generate_primes(fb_bound, primes);
    if (fb_size == 0)
    {
        fprintf(stderr, "Error: factor base generation failed\n");
        return 0;
    }
    
    relation_count = 0;
    matrix_rows = 0;
    int col_words = (fb_size + 63) / 64; // only algebraic side counted in parity
    int combo_words = (MAX_REL + 63) / 64;
    
    u128 m = int_root(n > 1 ? n - 1 : n, degree); // approximate
    
    uint64_t dep_mask[(MAX_REL + 63) / 64];
    int target_rel = fb_size + 16; // small overshoot to force a dependency sooner
    
    for (int k = 1; k <= window; k++)
    {
        if (relation_count >= MAX_REL || relation_count >= target_rel)
            break;
        
        u128 a = m + (u128)k;
        u128 algebraic = pow_u128(a, degree) + 1; // f(a) = a^d + 1
        
        Relation rel;
        memset(&rel, 0, sizeof(rel));
        rel.a_offset = k;
        
        // Rational side fixed to 1 (all exponents 0)
        memset(rel.r_exp, 0, sizeof(rel.r_exp));
        if (!factor_with_fb(algebraic, primes, &fb_size, rel.a_exp))
            continue;
        
        // Build row parity bits: algebraic columns [0, fb_size)
        uint64_t row[col_words];
        memset(row, 0, sizeof(row));
        for (int i = 0; i < fb_size; i++)
        {
            if (rel.a_exp[i] % 2 == 1)
            {
                row[i / 64] |= (uint64_t)1 << (i % 64);
            }
        }
        
        // Save relation
        relations[relation_count] = rel;
        
        // Build combo bits (identity row)
        uint64_t combo[combo_words];
        memset(combo, 0, sizeof(combo));
        combo[relation_count / 64] |= (uint64_t)1 << (relation_count % 64);
        
        if (insert_row(row, combo, col_words, combo_words, dep_mask))
        {
            u128 factor = attempt_dependency(dep_mask, combo_words, primes, fb_size, n);
            if (factor > 1 && factor < n)
                return factor;
        }
        relation_count++;
    }
    
    return 0;
}

// ============ CLI / demo ============

void run_demo()
{
    const char *demo_n_str = "815730722"; // 13^8 + 1 (small, finishes fast)
    u128 n = parse_u128(demo_n_str);
    int degree = 8;
    int fb = 200;
    int K = 5000; // search bound for k
    
    printf("SNFS Demo (toy) on n = ");
    print_u128(n);
    printf(" (degree=%d, B=%d, K=%d)\n\n", degree, fb, K);
    
    clock_t start = clock();
    u128 p = snfs_factor(n, degree, fb, K);
    clock_t mid = clock();
    double elapsed = (double)(mid - start) / CLOCKS_PER_SEC;
    
    if (p == 0 || p == n)
    {
        printf("SNFS toy failed, trying Pollard rho fallback...\n");
        p = pollard_rho_u128(n);
    }
    clock_t end = clock();
    double elapsed_total = (double)(end - start) / CLOCKS_PER_SEC;
    
    if (p == 0 || p == n)
    {
        printf("Failed to factor.\n");
        return;
    }
    u128 q = n / p;
    
    printf("Factors:\n  p = ");
    print_u128(p);
    printf("\n  q = ");
    print_u128(q);
    printf("\nSNFS time: %.4fs, total with fallback: %.4fs\n", elapsed, elapsed_total);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Usage: %s <n> [e] [degree] [B] [K]\n", argv[0]);
        printf("       %s --demo\n", argv[0]);
        return 1;
    }
    
    if (strcmp(argv[1], "--demo") == 0)
    {
        run_demo();
        return 0;
    }
    
    u128 n = parse_u128(argv[1]);
    u128 e = (argc >= 3) ? parse_u128(argv[2]) : 3;
    int degree = (argc >= 4) ? atoi(argv[3]) : 8;
    int fb = (argc >= 5) ? atoi(argv[4]) : 200;
    int K = (argc >= 6) ? atoi(argv[5]) : 5000; // k bound
    
    if (degree < 3 || degree > 12)
    {
        fprintf(stderr, "Degree must be between 3 and 12 for this toy.\n");
        return 1;
    }
    
    printf("SNFS (toy) Factorization\n");
    printf("n = ");
    print_u128(n);
    printf("\ne = ");
    print_u128(e);
    printf("\ndegree = %d, B = %d, K = %d\n\n", degree, fb, K);
    
    clock_t start = clock();
    u128 p = snfs_factor(n, degree, fb, K);
    clock_t mid = clock();
    double elapsed = (double)(mid - start) / CLOCKS_PER_SEC;
    
    if (p == 0 || p == n)
    {
        printf("SNFS toy failed, trying Pollard rho fallback...\n");
        p = pollard_rho_u128(n);
    }
    clock_t end = clock();
    double elapsed_total = (double)(end - start) / CLOCKS_PER_SEC;
    
    if (p == 0 || p == n)
    {
        printf("Failed to factor (try increasing B or K).\n");
        return 1;
    }
    
    u128 q = n / p;
    
    printf("Factors found:\n  p = ");
    print_u128(p);
    printf("\n  q = ");
    print_u128(q);
    printf("\nSNFS time: %.4fs, total with fallback: %.4fs\n\n", elapsed, elapsed_total);
    
    // Compute private key info if possible
    u128 phi = (p - 1) * (q - 1);
    if (gcd_u128(e, phi) == 1)
    {
        u128 d = mod_inverse_u128(e, phi);
        printf("phi(n) = ");
        print_u128(phi);
        printf("\nprivate exponent d = ");
        print_u128(d);
        printf("\n");
    }
    else
    {
        printf("e not coprime to phi(n), skipping d.\n");
    }
    
    return 0;
}
