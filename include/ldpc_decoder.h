#ifndef LDPC_DECODER_H
#define LDPC_DECODER_H

/**
 * @file ldpc_decoder.h
 * @brief LDPC Sum-Product (SPA) decoder and bit-wise LLR utilities.
 *
 * This header declares:
 *   - A standard LDPC decoder based on the Sum-Product Algorithm (SPA)
 *     in the log-likelihood ratio (LLR) domain
 *   - A helper function to convert symbol-wise likelihoods into
 *     bit-wise LLR values (for arbitrary modulation order E)
 *
 * The decoding algorithm follows the Tanner-graph message passing
 * schedule (flooding: CN update → VN update).
 *
 * Assumptions:
 *   - All LDPC operations use GF(2) arithmetic on the code side.
 *   - Code is systematic: codeword = [parity bits | information bits].
 *   - Channel LLRs follow the sign convention:
 *        LLR = log( P(y|x=+1) / P(y|x=-1) )
 */

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 *  LDPC Decoder: Sum-Product Algorithm (SPA)
 * ============================================================================
 *
 *  Performs iterative belief propagation on the Tanner graph defined by H.
 *
 *  Parameters:
 *      LLR[N]    : Input channel LLRs for each code bit
 *      ecc[N]    : Output decoded codeword bits (0/1)
 *      inf[K]    : Output decoded information bits (0/1)
 *      H[M][N]   : Parity-check matrix (entries are 0/1)
 *      M         : Number of parity-check equations
 *      N         : Codeword length
 *      K         : Information word length (systematic tail)
 *      max_iter  : Maximum number of SPA iterations
 *
 *  Notes:
 *      - Flooding schedule: check-node update → variable-node update.
 *      - Early termination if all parity checks are satisfied.
 *      - Assumes systematic code: info bits are extracted from
 *            ecc[N-K ... N-1]
 */
void ldpc_decode_spa(double *LLR, int *ecc, int *inf, int **H, int M, int N,
                     int K, int max_iter);

/* ============================================================================
 *  Compute bit-wise LLR from symbol-wise likelihoods
 * ============================================================================
 *
 *  Converts P(y_i | x = symbol_k) into bit-wise LLRs:
 *
 *      LLR[b + i*log2(E)] =
 *          log( Σ_{k : bit_b(k)=1} pyx[k][i] /
 *               Σ_{k : bit_b(k)=0} pyx[k][i] )
 *
 *  Parameters:
 *      pyx[E][N] : likelihoods for each symbol (k) and time index (i)
 *      E         : modulation order (must be a power of two)
 *      N         : number of symbols in the frame
 *      LLR       : output array of size N * log2(E)
 *
 *  Notes:
 *      - Supports any M-ary modulation (BPSK, QPSK, 8PSK, 16QAM, ...)
 *      - Symbol bit-labels are derived from binary representation of k
 */
void compute_llr_from_pyx(double **pyx, int E, int N, double *LLR);

#ifdef __cplusplus
}
#endif

#endif /* LDPC_DECODER_H */
