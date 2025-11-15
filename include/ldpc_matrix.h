/**
 * @file ldpc_matrix.h
 * @brief LDPC parity-check (H) and generator (G) matrix utilities.
 *
 * This header declares helper functions for LDPC code construction and
 * analysis:
 *
 *   1) Regular Gallager-type LDPC parity-check matrix generation (H)
 *   2) Systematic generator matrix construction (G) from H via Gaussian
 * elimination 3) 4-cycle counting for structural evaluation of LDPC codes
 *
 * All matrix operations are performed over GF(2), i.e., addition is XOR.
 */

#ifndef LDPC_MATRIX_H
#define LDPC_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* 1. LDPC Parity-Check Matrix Generation (Gallager Construction)             */
/* ========================================================================== */
/**
 * @brief Generate a regular (w_c, w_r) LDPC parity-check matrix H.
 *
 * Constructs an M×N parity-check matrix with:
 *    - Column weight  w_c  (each column has w_c ones)
 *    - Row weight     w_r  (each row has w_r ones)
 *
 * where:
 *    M = N * w_c / w_r
 *
 * The matrix is built using Gallager's block construction:
 *
 *   - H is divided into w_c row-blocks of size (M / w_c) × N
 *   - The first block has deterministically placed ones
 *   - Remaining (w_c − 1) blocks are random column permutations
 *     of the first block, preserving regularity
 *
 * @param H  Output parity-check matrix (allocated externally), size M × N
 * @param N  Codeword length (number of columns)
 * @param wc Column weight (number of ones per column)
 * @param wr Row weight (number of ones per row)
 */
void generate_Hmatrix(int **H, int N, int wc, int wr);

/* ========================================================================== */
/* 2. Systematic Generator Matrix Construction                                */
/* ========================================================================== */
/**
 * @brief Construct a systematic generator matrix G from an LDPC H matrix.
 *
 * Given a parity-check matrix H of size M×N, this routine derives a
 * systematic generator matrix G of size K×N, where:
 *
 *   M = N * w_c / w_r
 *   K = N − M
 *
 * Method (all operations in GF(2)):
 *
 *   1) Form the augmented matrix:
 *          X = [ H^T | I_N ]      (size N × (M+N))
 *
 *   2) Perform Gaussian elimination on X to obtain a systematic form
 *      for the left block (H^T part).
 *
 *   3) Continue elimination on the right block while allowing column
 *      swaps. Column swaps that affect X are also mirrored into H so
 *      that the constraint H·G^T = 0 is preserved.
 *
 *   4) Extract the generator matrix G from the lower K rows of the
 *      right block:
 *          G = X[M..N-1][M..M+N-1]
 *
 * The resulting G is systematic (typically G = [P | I_K] up to
 * column permutations).
 *
 * @param H  Parity-check matrix (M × N). May be modified in-place due
 *           to column permutations during elimination.
 * @param G  Output generator matrix (allocated externally), size K × N
 * @param N  Codeword length
 * @param wc Column weight used for H
 * @param wr Row weight used for H
 */
void generate_Gmatrix(int **H, int **G, int N, int wc, int wr);

/* ========================================================================== */
/* 3. Structural Analysis: 4-Cycle Counting                                   */
/* ========================================================================== */
/**
 * @brief Count the number of 4-cycles in an LDPC parity-check matrix H.
 *
 * A 4-cycle arises when two variable nodes (columns) share two or more
 * common check nodes (rows). Short cycles (especially length-4) degrade
 * the performance of message-passing decoders such as SPA/BP.
 *
 * Algorithm (conceptual):
 *   - For each column, record the indices of the rows where H(row, col) = 1
 *   - For every unordered pair of columns, count how many rows they share
 *   - If two columns share `shared >= 2` check nodes, they contribute
 *     nC2 = shared! / (2!(shared − 2)!) distinct 4-cycles
 *
 * This function is intended for structural analysis only; it does not
 * modify H.
 *
 * @param H   Parity-check matrix, size M × N (M = N * w_c / w_r)
 * @param N   Codeword length
 * @param wc  Column weight
 * @param wr  Row weight
 *
 * @return    Total number of 4-cycles detected in H.
 */
int count_floop(int **H, int N, int wc, int wr);

#ifdef __cplusplus
}
#endif

#endif /* LDPC_MATRIX_H */
