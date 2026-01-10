# RSA Interactive

Minimal C implementation that demonstrates RSA key generation, encryption, and CRT-based decryption, printing each step.

## Files
- README.md: technical notes.
- rsa_interactive.c: full source for keygen/encrypt/decrypt.
- trial_division.c / pollards_rho.c: basic factorization demos.
- snfs.c: toy Special NFS-style factorer with fallback to Pollard rho.

## Requirements
- gcc (or any C11 compiler).
- Basic POSIX libs (stdio.h, stdlib.h).

## Build and run
```bash
gcc rsa_interactive.c -o rsa_interactive
./rsa_interactive

gcc trial_division.c -o trial_division
gcc pollards_rho.c -o pollards_rho
gcc snfs.c -o snfs
```
The binary asks for a message (up to 1023 chars), encrypts per character, then decrypts with CRT and compares to the original.

### Factorization demos
- Trial division: `./trial_division <n>`
- Pollard’s rho: `./pollards_rho <n>`
- Toy SNFS (special-form n): `./snfs <n> [e] [degree] [B] [K]`
  - Example (works fast): `./snfs 815730722 3 8 200 5000` (`n = 13^8 + 1`)
  - For larger special forms (e.g., `614^8 + 1 = 20199795332516287488257`), the toy SNFS is unlikely to finish; you’ll need a real NFS implementation (msieve, cado-nfs) or accept a Pollard fallback.

## Program flow
1. Uses fixed exponent `e = 3`.
2. Picks pseudo-random 16-bit primes `p` and `q` via `rand()` and naive `ifprime`, ensuring `gcd(e, phi) = 1`.
3. Computes `n = p * q`, `phi = (p-1)*(q-1)`, and modular inverse `d = e^-1 mod phi`.
4. Reads user message.
5. Encrypts each byte: `c = m^e mod n` (`modpow_encrypt`).
6. Decrypts each `c` with CRT: `dP = d mod (p-1)`, `dQ = d mod (q-1)`, `qInv = q^-1 mod p`, then reconstructs `m`.
7. Compares original vs decrypted and reports the result.

## Key functions
- `setprimes`: generates valid `p`, `q`, `n`, `phi` for chosen `e`.
- `findD`: computes `d`, modular inverse of `e` mod `phi`.
- `encrypt_text`: per-byte encryption via `modpow_encrypt`, tracking length.
- `decrypt_text`: CRT-based decryption using `modpow_decrypt`, `inverse`, and `m = m2 + h * q`.
- `gcd`, `ifprime`, `getprime`: utilities for primality and coprimality.

## Limitations (educational only)
- Very small keys (16-bit); no real cryptographic security.
- Primes from `rand()` and naive primality test; predictable and repeatable.
- `modpow_encrypt`/`modpow_decrypt` use linear-time exponentiation; slow for larger exponents.
- No overflow protection; `unsigned long long` may not hold larger powers.
- ASCII-only; no multibyte handling.
- Minimal I/O error handling.

## Quick example
```
$ ./rsa_interactive
========================================
   Interactive RSA Encryption System
========================================

Step 1: Generating RSA Keys...
...
Enter your message: hello rsa
...
✓ SUCCESS: Decryption matches original!
```
