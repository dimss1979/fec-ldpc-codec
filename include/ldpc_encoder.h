/**
 * @file ldpc_encoder.h
 * @brief LDPC encoder using generator-matrix multiplication over GF(2).
 *
 * This module provides a simple and general-purpose LDPC encoder based on:
 *
 *      c = u · G      (mod 2)
 *
 * where:
 *    - u : K-bit information vector
 *    - G : K×N generator matrix (systematic form preferred)
 *    - c : N-bit encoded codeword
 *
 * All arithmetic is over GF(2):
 *    - addition → XOR
 *    - multiplication → AND
 */

#ifndef LDPC_ENCODER_H
#define LDPC_ENCODER_H

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 *  Function: ldpc_encode
 *  --------------------------------------------------------------------------
 *  Encode information bits using a generator matrix G (GF(2)).
 *
 *  Parameters:
 *      ecc : Output codeword (int array of length N)
 *      inf : Input information bits (int array of length K, values 0/1)
 *      G   : Generator matrix G[j][i] (size K × N)
 *      N   : Codeword length
 *      K   : Information length
 *
 *  Operation:
 *      ecc[i] = XOR_{j=0..K-1} ( inf[j] & G[j][i] )
 *
 *  Notes:
 *      - Assumes G is already systematic (e.g., [P | I]) but not required.
 *      - All arithmetic is performed in GF(2) via XOR and AND.
 * ============================================================================
 */
void ldpc_encode(int *ecc, const int *inf, int **G, int N, int K);

#ifdef __cplusplus
}
#endif

#endif /* LDPC_ENCODER_H */
