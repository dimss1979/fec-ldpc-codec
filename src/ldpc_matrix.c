/**
 * @file ldpc_matrix.c
 * @brief LDPC parity-check matrix (H) and generator matrix (G) utilities
 *
 * This module provides:
 *   - Regular (w_c, w_r) Gallager-type parity-check matrix generation
 *   - Systematic generator matrix construction via Gaussian elimination
 *   - 4-cycle counting for structural LDPC code evaluation
 *
 * All operations are over GF(2): addition = XOR, multiplication = AND.
 *
 * The implementation follows a classical LDPC workflow:
 *
 *    1) Construct H (parity-check matrix)
 *    2) Derive G (systematic generator matrix)
 *    3) Analyze H (cycle-4 count, degree constraints, etc.)
 *
 * This file supports LDPC research, coding-theory education, and BER
 * simulation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ldpc_matrix.h"

/* ========================================================================== */
/* 1. Gallager Regular LDPC Parity-Check Matrix Generation                    */
/* -------------------------------------------------------------------------- */
/*
 * Construct an M×N regular LDPC matrix (column-weight wc, row-weight wr).
 *
 *   M = N * (wc / wr)
 *
 * Gallager's method:
 *   - H is divided into wc row-blocks
 *   - The first block is constructed deterministically
 *   - Remaining blocks are random column-permutations of the first block
 *
 * Example (wc = 3, wr = 6):
 *       [ B0 ]
 *   H = [ P1(B0) ]
 *       [ P2(B0) ]      where Pk is a different random permutation
 *
 * This yields a sparse, regular LDPC code with controllable 4-cycle behavior.
 */
/* ========================================================================== */
void generate_Hmatrix(int **H, int N, int wc, int wr) {
  int i, j, k;
  int M = (N * wc) / wr;   /* number of check equations */
  int block_rows = M / wc; /* rows per Gallager block   */

  int *perm = (int *)malloc(N * sizeof(int));
  if (!perm) {
    fprintf(stderr, "malloc failed in generate_Hmatrix\n");
    exit(1);
  }

  /* initialize H to zero */
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++)
      H[i][j] = 0;

  /* ---------------- Block 0: deterministic spacing of '1's ---------------- */
  for (i = 0; i < block_rows; i++)
    for (j = i * wr; j < (i + 1) * wr; j++)
      H[i][j] = 1;

  /* ----------- Block 1 ... wc-1: apply random column permutations ----------
   */
  for (i = 1; i < wc; i++) {

    /* identity permutation */
    for (j = 0; j < N; j++)
      perm[j] = j;

    /* Fisher-Yates shuffle */
    for (j = 0; j < N; j++) {
      int r = rand() % N;
      int tmp = perm[j];
      perm[j] = perm[r];
      perm[r] = tmp;
    }

    /* apply permuted columns into block i */
    for (j = block_rows * i; j < block_rows * (i + 1); j++)
      for (k = 0; k < N; k++)
        H[j][k] = H[j - block_rows * i][perm[k]];
  }

  free(perm);
}

/* ========================================================================== */
/* 2. Systematic Generator Matrix Construction (G from H)                     */
/* -------------------------------------------------------------------------- */
/*
 * Construct systematic generator matrix G (size K×N) from parity-check matrix
 * H.
 *
 * Definitions:
 *      M = number of parity-check equations
 *      N = code length
 *      K = N − M  (information bits)
 *
 * Method:
 *    Form augmented matrix:
 *
 *        X = [ H^T | I_N ]     (size N × (M+N))
 *
 *    Perform Gaussian elimination on X to obtain systematic form:
 *
 *        X → [ I_M | A ]
 *             [  0  | G ]
 *
 *    The lower block yields the generator matrix G:
 *
 *        G = X[M..N-1][M..M+N-1]
 *
 * Notes:
 *   - All operations are in GF(2)
 *   - Some elimination steps require column swaps
 *   - Swapping columns in X corresponds to swapping columns in H
 *     to maintain the H·G^T = 0 constraint
 */
/* ========================================================================== */
void generate_Gmatrix(int **H, int **G, int N, int wc, int wr) {
  int M = (N * wc) / wr;

  int i, j, k, l;

  /* X is the augmented matrix: [H^T | I] */
  int **X = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++)
    X[i] = (int *)malloc((M + N) * sizeof(int));

  int *row_buf = (int *)malloc((M + N) * sizeof(int));
  int *col_buf = (int *)malloc(N * sizeof(int));
  int *Hcol_buf = (int *)malloc(M * sizeof(int));

  if (!X || !row_buf || !col_buf || !Hcol_buf) {
    fprintf(stderr, "malloc failed in generate_Gmatrix\n");
    exit(1);
  }

  /* --------------------- Step 1: Build [H^T | I] ------------------------ */
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++)
      X[i][j] = H[j][i]; /* left block */
    for (j = M; j < M + N; j++)
      X[i][j] = (i == (j - M)) ? 1 : 0; /* right block = I */
  }

  /* -------- Step 2: Gaussian elimination on left block (H^T part only) --- */
  for (j = 0; j < M; j++) {

    /* If pivot missing → row swap OR column swap within X */
    if (X[j][j] == 0) {

      int pivot_found = 0;
      for (i = j + 1; i < N; i++) {
        if (X[i][j] == 1) {
          memcpy(row_buf, X[i], (M + N) * sizeof(int));
          memcpy(X[i], X[j], (M + N) * sizeof(int));
          memcpy(X[j], row_buf, (M + N) * sizeof(int));
          pivot_found = 1;
          break;
        }
      }

      /* If pivot still not found → swap columns inside X */
      if (!pivot_found) {
        for (k = M + N - 1; k > j; k--) {
          if (X[j][k] == 1) {
            for (l = 0; l < N; l++) {
              col_buf[l] = X[l][k];
              X[l][k] = X[l][j];
              X[l][j] = col_buf[l];
            }
            break;
          }
        }
      }
    }

    /* Row elimination (GF(2)) */
    for (i = 0; i < N; i++)
      if (i != j && X[i][j] == 1)
        for (k = 0; k < M + N; k++)
          X[i][k] ^= X[j][k];
  }

  /* ------------- Step 3: Elimination on right block with H updates ------- */
  for (j = 2 * M; j < M + N; j++) {

    int pivot_row = j - M;

    if (X[pivot_row][j] == 0) {

      int found = 0;
      for (i = pivot_row + 1; i < N; i++) {
        if (X[i][j] == 1) {
          memcpy(row_buf, X[i], (M + N) * sizeof(int));
          memcpy(X[i], X[pivot_row], (M + N) * sizeof(int));
          memcpy(X[pivot_row], row_buf, (M + N) * sizeof(int));
          found = 1;
          break;
        }
      }

      /* Still missing pivot → swap columns → update H accordingly */
      if (!found) {
        for (k = M + N - 1; k > M - 1; k--) {
          if (X[pivot_row][k] == 1) {

            /* swap columns in X */
            for (l = 0; l < N; l++) {
              col_buf[l] = X[l][k];
              X[l][k] = X[l][j];
              X[l][j] = col_buf[l];
            }

            /* mirror swap inside H (only H columns affected) */
            for (l = 0; l < M; l++) {
              Hcol_buf[l] = H[l][k - M];
              H[l][k - M] = H[l][j - M];
              H[l][j - M] = Hcol_buf[l];
            }

            break;
          }
        }
      }
    }

    /* eliminate other rows */
    for (i = 0; i < N; i++)
      if (i != pivot_row && X[i][j] == 1)
        for (k = 0; k < M + N; k++)
          X[i][k] ^= X[pivot_row][k];
  }

  /* ----------------------- Step 4: Extract G (K×N) ----------------------- */
  for (i = M; i < N; i++)
    for (j = M; j < M + N; j++)
      G[i - M][j - M] = X[i][j];

  /* cleanup */
  for (i = 0; i < N; i++)
    free(X[i]);
  free(X);
  free(row_buf);
  free(col_buf);
  free(Hcol_buf);
}

/* ========================================================================== */
/* 3. 4-Cycle Counting in LDPC Parity-Check Matrix */
/* -------------------------------------------------------------------------- */
/*
 * A 4-cycle exists when two variable nodes share ≥2 check nodes.
 * Short cycles harm message-passing performance (SPA/BP).
 *
 * This function:
 *   - Builds adjacency lists of check nodes for each variable
 *   - Counts the number of pairs of columns sharing ≥2 rows
 *   - Counts multiplicities using nC2 = shared! / (2!(shared-2)!)
 *
 * A good LDPC code strives to minimize such 4-cycles.
 */
/* ========================================================================== */
static int factorial(int n) {
  int f = 1;
  for (int i = 1; i <= n; i++)
    f *= i;
  return f;
}

int count_floop(int **H, int N, int wc, int wr) {
  int M = (N * wc) / wr;

  /* store row positions of '1's for each column */
  int **var_nodes = (int **)malloc(N * sizeof(int *));
  for (int i = 0; i < N; i++)
    var_nodes[i] = (int *)malloc(wc * sizeof(int));

  /* build adjacency list: var_nodes[j][*] = rows where H[row][j] = 1 */
  for (int j = 0; j < N; j++) {
    int idx = 0;
    for (int i = 0; i < M; i++) {
      if (H[i][j])
        var_nodes[j][idx++] = i;
      if (idx == wc)
        break;
    }
  }

  int floop = 0;

  /* Count column pairs sharing ≥2 check nodes (4-cycles) */
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {

      int shared = 0;
      for (int k = 0; k < wc; k++)
        for (int l = 0; l < wc; l++)
          if (var_nodes[i][k] == var_nodes[j][l])
            shared++;

      if (shared >= 2)
        floop += factorial(shared) / (2 * factorial(shared - 2));
    }
  }

  /* cleanup */
  for (int i = 0; i < N; i++)
    free(var_nodes[i]);
  free(var_nodes);

  return floop;
}
