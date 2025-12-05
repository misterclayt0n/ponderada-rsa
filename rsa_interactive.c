#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>

#define MAX_VALUE 65535
#define E_VALUE 3
#define MAX_TEXT_LENGTH 1024

uint32_t findD(uint16_t e, uint32_t phi)
{
	uint32_t eprev, dprev, d = 1, etemp, dtemp;
	
	eprev = phi, dprev = phi;
	while (e != 1)
	{
		etemp = e;
		dtemp = d;
		e = eprev - eprev / etemp * e;
		d = dprev - eprev / etemp * d;
		eprev = etemp;
		dprev = dtemp;
		while (d < 0)
			d += phi;
	}
	return d;
}

int ifprime(uint16_t n)
{
	for (uint16_t i = 2; i <= n / 2; i++)
	{
		if (n % i == 0)
			return 0;
	}
	return 1;
}

uint16_t gcd(uint16_t num1, uint32_t num2)
{
	uint16_t temp;
	if (num1 > num2)
	{
		temp = num1;
		num1 = num2;
		num2 = temp;
	}
	for (uint16_t i = num1; i > 0; i--)
	{
		if (num1 % i == 0 && num2 % i == 0)
			return i;
	}
	return 1;
}

uint16_t getprime()
{
	uint16_t n;
	do
	{
		n = rand() % MAX_VALUE + 5;
	} while (!ifprime(n));
	return n;
}

void setprimes(uint16_t e, uint16_t *p, uint16_t *q, uint32_t *n, uint32_t *phi)
{
	do
	{
		*p = getprime();
		do
			*q = getprime();
		while(*p == *q);
		
		*n = *p * *q;
		*phi = *n - *p - *q + 1;
	} while (gcd(e, *phi) != 1);
}

unsigned long long int modpow_encrypt(int base, int power, int mod)
{
	unsigned long long int result = 1;
	for (int i = 0; i < power; i++)
		result = (result * base) % mod;
	return result;
}

void encrypt_text(const char *plaintext, unsigned long long int *ciphertext, int *cipher_len, uint32_t n, uint16_t e)
{
	*cipher_len = 0;
	for (int i = 0; plaintext[i] != '\0' && plaintext[i] != '\n'; i++)
	{
		ciphertext[i] = modpow_encrypt((unsigned char)plaintext[i], e, n);
		(*cipher_len)++;
	}
}

unsigned long long int modpow_decrypt(unsigned long long int base, int power, int mod)
{
	unsigned long long int result = 1;
	for (int i = 0; i < power; i++)
		result = (result * base) % mod;
	return result;
}

int inverse(int a, int mod)
{
	int aprev, iprev, i = 1, atemp, itemp;
	
	aprev = mod, iprev = mod;
	while (a != 1)
	{
		atemp = a;
		itemp = i;
		a = aprev - aprev / atemp * a;
		i = iprev - aprev / atemp * i;
		aprev = atemp;
		iprev = itemp;
		while (i < 0)
			i += mod;
	}
	return i;
}

void decrypt_text(unsigned long long int *ciphertext, int cipher_len, char *plaintext, 
                  uint32_t n, uint32_t d, uint16_t p, uint16_t q)
{
	unsigned long long int dP, dQ, m1, m2;
	int qInv, m1m2, h, m;
	
	dP = d % (p - 1);
	dQ = d % (q - 1);
	qInv = inverse(q, p);
	
	for (int i = 0; i < cipher_len; i++)
	{
		m1 = modpow_decrypt(ciphertext[i], dP, p);
		m2 = modpow_decrypt(ciphertext[i], dQ, q);
		m1m2 = m1 - m2;
		if (m1m2 < 0)
			m1m2 += p;
		h = (qInv * m1m2) % p;
		m = m2 + h * q;
		plaintext[i] = (char)m;
	}
	plaintext[cipher_len] = '\0';
}

int main()
{
	uint16_t e = E_VALUE, p, q;
	uint32_t n, phi, d;
	char plaintext[MAX_TEXT_LENGTH];
	unsigned long long int ciphertext[MAX_TEXT_LENGTH];
	char decrypted[MAX_TEXT_LENGTH];
	int cipher_len;
	
	srand(time(NULL));
	
	printf("RSA Encryption System\n\n");
	
	setprimes(e, &p, &q, &n, &phi);
	d = findD(e, phi);
	
	printf("Keys generated:\n");
	printf("  p = %"PRIu16", q = %"PRIu16"\n", p, q);
	printf("  n = %"PRIu32", phi = %"PRIu32"\n", n, phi);
	printf("  e = %"PRIu16", d = %"PRIu32"\n\n", e, d);
	
	printf("Enter message: ");
	if (fgets(plaintext, MAX_TEXT_LENGTH, stdin) == NULL)
	{
		printf("Error reading input\n");
		return 1;
	}
	
	size_t len = strlen(plaintext);
	if (len > 0 && plaintext[len-1] == '\n')
		plaintext[len-1] = '\0';
	
	encrypt_text(plaintext, ciphertext, &cipher_len, n, e);
	
	printf("\nCiphertext: ");
	for (int i = 0; i < cipher_len; i++)
		printf("%llu ", ciphertext[i]);
	printf("\n");
	
	decrypt_text(ciphertext, cipher_len, decrypted, n, d, p, q);
	
	printf("\nOriginal:  \"%s\"\n", plaintext);
	printf("Decrypted: \"%s\"\n", decrypted);
	printf("Status: %s\n", strcmp(plaintext, decrypted) == 0 ? "OK" : "FAILED");
	
	return 0;
}
