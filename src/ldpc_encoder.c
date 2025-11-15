/**
 * @file ldpc_encoder.c
 * @brief LDPC encoder using systematic generator matrix (G).
 *
 * This module implements LDPC encoding via matrix multiplication in GF(2):
 *
 *      c = u · G      (mod 2)
 *
 * where:
 *    - u   : K-bit information vector
 *    - G   : K×N systematic generator matrix (typically [P | I])
 *    - c   : N-bit encoded codeword
 *
 * All operations are XOR-based (GF(2)).
 */

#include "ldpc_encoder.h"

/* ========================================================================
 *  LDPC Encoder (Generator-Matrix Based)
 * ------------------------------------------------------------------------
 *  ecc[i] = XOR_{j=0..K-1} ( inf[j] AND G[j][i] )
 *
 *  Parameters:
 *      ecc : Output codeword (length N)
 *      inf : Input information bits (length K)
 *      G   : Generator matrix, size K × N (row-major)
 *      N   : Codeword length
 *      K   : Information length
 *
 *  Notes:
 *      - Assumes G is already in systematic form (e.g., [P | I]).
 *      - Encoding complexity: O(KN) bit operations.
 *      - All arithmetic is in GF(2): additions → XOR, multiplications → AND.
 * ======================================================================== */
void ldpc_encode(int *ecc, const int *inf, int **G, int N, int K) {
  for (int i = 0; i < N; i++) {

    int acc = 0; /* accumulator for column i */

    /* Multiply information vector by column i of G */
    for (int j = 0; j < K; j++) {
      acc ^= (inf[j] & G[j][i]);
    }

    ecc[i] = acc;
  }
}
