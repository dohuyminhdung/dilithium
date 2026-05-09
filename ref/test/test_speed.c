#include <stdint.h>
#include "../sign.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../params.h"
// #include "cpucycles.h"
#include "speed_print.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../randombytes.h"
#include "../fips202.h"
#include "../packing.h"

#define NTESTS 10000
#define SO_FILE 10
#define ExpandA       0
#define ExpandS       1
#define ExpandMask    2
#define NTT           3
#define INTT          4
#define MultiplyPoly  5
#define SampleInBall  6
#define Keygen        7
#define Sign          8
#define Verify        9

static keccak_state rngstate = {{0x1F, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1ULL << 63), 0, 0, 0, 0}, SHAKE128_RATE};
void randombytes_benchmark(uint8_t *x,size_t xlen) {
  shake128_squeeze(x, xlen, &rngstate);
}

void poly_write(FILE *fp, const poly *p)
{
    fprintf(fp, "[");

    for(int i = 0; i < N; i++) {
        fprintf(fp, "%d", p->coeffs[i]);

        if (i != N - 1) {
            fprintf(fp, ", ");
        }
    }

    fprintf(fp, "]\n");
}

int main(void)
{
  unsigned int i;
  size_t siglen;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sig[CRYPTO_BYTES];
  uint8_t seed[CRHBYTES];
  uint8_t msg[256];

  polyvecl mat[K];
  poly *a = &mat[0].vec[0];
  poly *b = &mat[0].vec[1];
  poly *c = &mat[0].vec[2];
  
  //benchmark in ARM
  struct timespec t1, t2;
  double elapsed_us = 0;
  double average_us = 0;
  double average_us1 = 0;
  double average_us2 = 0;

  //dump input
  FILE *files[SO_FILE];
  const char *tenFiles[SO_FILE] = {
        "ExpandA.txt",
        "ExpandS.txt",
        "ExpandMask.txt",
        "NTT.txt",
        "INTT.txt",
        "MultiplyPoly.txt",
        "SampleInBall.txt",
        "Keygen.txt",
        "Sign.txt",
        "Verify.txt"
    };

    for(i = 0; i < SO_FILE; i++) {
        files[i] = fopen(tenFiles[i], "w");
        if (files[i] == NULL) {
            printf("Khong the mo file %s\n", tenFiles[i]);
            for(unsigned int j = 0; j < i; j++) {
                fclose(files[j]);
            }
            return 1;
        }
    }


  for(i = 0; i < NTESTS; ++i) {
    randombytes_benchmark(seed, 32);
    for(int j = 0; j < 32; j++) {
      fprintf(files[ExpandA], "%02X", seed[j]);
    }
    fprintf(files[ExpandA], "\n");
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    polyvec_matrix_expand(mat, seed);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
  }
  average_us = average_us / NTESTS;
  printf("polyvec_matrix_expand: %f\n", average_us);
  fclose(files[ExpandA]);
  
  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    randombytes_benchmark(seed, 64);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    poly_uniform_eta(a, seed, 0);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
    for(int j = 0; j < 64; j++) {
      fprintf(files[ExpandS], "%02X", seed[j]);
    }
    fprintf(files[ExpandS], "\n");
  }
  average_us = average_us / NTESTS * 15;
  printf("poly_uniform_eta: %f\n", average_us);
  fclose(files[ExpandS]);

  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    randombytes_benchmark(seed, 64);
        for(int j = 0; j < 64; j++) {
      fprintf(files[ExpandMask], "%02X", seed[j]);
    }
    fprintf(files[ExpandMask], "\n");
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    poly_uniform_gamma1(a, seed, 0);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
  }
  average_us = average_us / NTESTS * 7;
  printf("poly_uniform_gamma1: %f\n", average_us);
  fclose(files[ExpandMask]);

  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    poly_ntt(a);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
    poly_write(files[NTT], a);
  }
  average_us = average_us / NTESTS;
  printf("poly_ntt: %f\n", average_us);
  fclose(files[NTT]);

  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    poly_invntt_tomont(a);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
    poly_write(files[INTT], a);
  }
  average_us = average_us / NTESTS;
  printf("poly_invntt_tomont: %f\n", average_us);
  fclose(files[INTT]);

  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    poly_pointwise_montgomery(c, a, b);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
    poly_write(files[MultiplyPoly], a);
    poly_write(files[MultiplyPoly], b);
    fprintf(files[MultiplyPoly], "\n");
  }
  average_us = average_us / NTESTS;
  printf("poly_pointwise_montgomery: %f\n", average_us);
  fclose(files[MultiplyPoly]);

  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    randombytes_benchmark(seed, 64);
    for(int j = 0; j < 64; j++) {
      fprintf(files[SampleInBall], "%02X", seed[j]);
    }
    fprintf(files[SampleInBall], "\n");
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    poly_challenge(c, seed);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
    elapsed_us = (t2.tv_sec - t1.tv_sec) * 1e6 +
                  (t2.tv_nsec - t1.tv_nsec) / 1e3;
    average_us = average_us + elapsed_us;
  }
  average_us = average_us / NTESTS;
  printf("poly_challenge: %f\n", average_us);
  fclose(files[SampleInBall]);

  average_us = 0;
  for(i = 0; i < NTESTS; ++i) {
    //KEY GEN
    fprintf(files[Keygen], "\n");
    crypto_sign_keypair_benchmark(pk, sk, seed, &elapsed_us);
    for(int j = 0; j < 32; j++) {
      fprintf(files[Keygen], "%02X", seed[j]);
    }
    average_us = average_us + elapsed_us;

    //get random msg
    randombytes_benchmark(seed, 64);

    //SIGN
    //sk
    fprintf(files[Sign], "Test %d:\n", i);
    fprintf(files[Sign], "secret key\n");
    for(int j = 0; j < CRYPTO_SECRETKEYBYTES; j++) {
      fprintf(files[Sign], "%02X", sk[j]);
    }
    fprintf(files[Sign], "\n");
    crypto_sign_signature_benchmark(sig, &siglen, seed, 64, NULL, 0, sk, msg, &elapsed_us);
    average_us1 = average_us1 + elapsed_us;
    //msg
    fprintf(files[Sign], "message\n");
    for(int j = 0; j < CRHBYTES + 2; j++) {
      fprintf(files[Sign], "%02X", msg[j]);
    }
    fprintf(files[Sign], "\n");

    //VERIFY
    fprintf(files[Verify], "Test %d:\n", i);
    //pk
    fprintf(files[Verify], "public key\n");
    for(int j = 0; j < CRYPTO_PUBLICKEYBYTES; j++) {
      fprintf(files[Verify], "%02X", pk[j]);
    }
    fprintf(files[Verify], "\n");
    //sig
    fprintf(files[Verify], "signature\n");
    for(int j = 0; j < CRYPTO_BYTES; j++) {
      fprintf(files[Verify], "%02X", sig[j]);
    }
    fprintf(files[Verify], "\n");
    crypto_sign_verify_benchmark(sig, CRYPTO_BYTES, seed, CRHBYTES, NULL, 0, pk, msg, &elapsed_us);
    //msg
    fprintf(files[Verify], "message\n");
    for(int j = 0; j < CRHBYTES + 2; j++) {
      fprintf(files[Verify], "%02X", msg[j]);
    }
    fprintf(files[Verify], "\n");
    average_us2 = average_us2 + elapsed_us;
  }
  average_us = average_us / NTESTS;
  printf("Keypair: %f\n", average_us);
  fclose(files[Keygen]);

  average_us1 = average_us1 / NTESTS;
  printf("Sign: %f\n", average_us1);
  fclose(files[Sign]);
  
  average_us2 = average_us2 / NTESTS;
  printf("Verify: %f \n", average_us2);
  fclose(files[Verify]);
  return 0;
}
