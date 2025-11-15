/**
 * @file ldpc_decoder.c
 * @brief LDPC Sum-Product (SPA) decoder and LLR computation utilities.
 *
 * This module provides:
 *   - A standard LDPC decoder based on the Sum-Product Algorithm (SPA)
 *     operating in the log-likelihood ratio (LLR) domain
 *   - A helper function to compute bit-wise LLRs from per-symbol likelihoods
 *
 * The SPA implementation uses:
 *   - Flooding schedule (check-node update → variable-node update)
 *   - LLR-domain message passing over the Tanner graph defined by H
 *   - Early stopping when all parity checks are satisfied (H·c^T = 0)
 *
 * Assumptions:
 *   - Code is systematic: codeword = [parity | info]
 *   - All LDPC operations are in GF(2) on the codeword side
 */

#include "ldpc_decoder.h"
#include <math.h>
#include <stdlib.h>

/* ========================================================================== */
/* Helper: sign(x)                                                            */
/* ========================================================================== */
/**
 * @brief Return the sign of x (+1 or -1), treating zero as positive.
 */
static inline double sign_val(double x) { return (x >= 0.0) ? 1.0 : -1.0; }

/* ========================================================================== */
/* Helper: spf(x) = log((e^x + 1)/(e^x - 1)), with safe clipping              */
/* ========================================================================== */
/**
 * @brief Check-node nonlinearity used in SPA (LLR-domain).
 *
 * For numerical stability, the input x is clipped into [1e-7, 30].
 * The function is used as an approximation to the standard SPA
 * check-node operation in the log-domain.
 */
static inline double spf(double x) {
  if (x < 1e-7)
    x = 1e-7;
  else if (x > 30.0)
    x = 30.0;

  return log((exp(x) + 1.0) / (exp(x) - 1.0));
}

/* ========================================================================== */
/* Sum-Product Algorithm (SPA) LDPC Decoder                                   */
/* ========================================================================== */
/**
 * @brief LDPC decoding using the Sum-Product Algorithm (LLR-domain, flooding).
 *
 * Tanner graph:
 *   - H: M×N parity-check matrix
 *   - Variable nodes: N
 *   - Check nodes   : M
 *
 * Message notation:
 *   - u[i][j] : message from variable node j → check node i (extrinsic LLR)
 *   - v[i][j] : message from check node i → variable node j (extrinsic LLR)
 *
 * Decoding steps per iteration:
 *   1) Check-node update:
 *        v[i][j] = f({ LLR[j'] + u[i][j'] | j' ≠ j, H[i][j']=1 })
 *   2) Variable-node update:
 *        u[i][j] = Σ_{i'≠i} v[i'][j]
 *   3) A-posteriori LLR:
 *        L_post[j] = LLR[j] + Σ_i v[i][j]
 *      → hard decision ecc[j] = (L_post[j] >= 0) ? 1 : 0
 *   4) Parity check:
 *        If H·ecc^T = 0, stop early.
 *
 * Finally, the information part is extracted assuming:
 *      codeword = [parity (N-K bits) | info (K bits)]
 *
 * @param LLR      Input channel LLRs for each code bit (length N)
 * @param ecc      Output decoded codeword bits (length N, 0/1)
 * @param inf      Output decoded information bits (length K, 0/1)
 * @param H        Parity-check matrix (M×N), entries in {0,1}
 * @param M        Number of parity-check equations (rows of H)
 * @param N        Codeword length (columns of H)
 * @param K        Information length (systematic part length)
 * @param max_iter Maximum number of SPA iterations
 */
void ldpc_decode_spa(double *LLR, int *ecc, int *inf, int **H, int M, int N,
                     int K, int max_iter) {
  int i, j, k, iter;

  /* ------------------------------------------------------------------ */
  /* Build adjacency lists for check and variable nodes                 */
  /* ------------------------------------------------------------------ */
  int **check_node =
      (int **)malloc(M * sizeof(int *)); /* neighbors of check i */
  for (i = 0; i < M; i++) {
    check_node[i] = (int *)malloc(N * sizeof(int));
  }

  int *deg_c = (int *)calloc(M, sizeof(int)); /* degree of each check node    */
  int *deg_v = (int *)calloc(N, sizeof(int)); /* degree of each variable node */

  /* Count degrees of check nodes */
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      if (H[i][j])
        deg_c[i]++;
    }
  }

  /* Count degrees of variable nodes */
  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      if (H[i][j])
        deg_v[j]++;
    }
  }

  /* Variable-node adjacency list: variable_node[j][*] = list of check indices
   */
  int **variable_node = (int **)malloc(N * sizeof(int *));
  for (j = 0; j < N; j++) {
    variable_node[j] = (int *)malloc(deg_v[j] * sizeof(int));
  }

  /* Fill check_node[i][*] = columns j for which H[i][j] = 1 */
  for (i = 0; i < M; i++) {
    int idx = 0;
    for (j = 0; j < N; j++) {
      if (H[i][j]) {
        check_node[i][idx++] = j;
      }
    }
  }

  /* Fill variable_node[j][*] = rows i for which H[i][j] = 1 */
  for (j = 0; j < N; j++) {
    int idx = 0;
    for (i = 0; i < M; i++) {
      if (H[i][j]) {
        variable_node[j][idx++] = i;
      }
    }
  }

  /* ------------------------------------------------------------------ */
  /* Allocate message arrays u[i][j] (V→C) and v[i][j] (C→V)            */
  /* For simplicity, these are allocated as M×N dense matrices.         */
  /* ------------------------------------------------------------------ */
  double **u = (double **)malloc(M * sizeof(double *));
  double **v = (double **)malloc(M * sizeof(double *));
  for (i = 0; i < M; i++) {
    u[i] = (double *)calloc(N, sizeof(double));
    v[i] = (double *)calloc(N, sizeof(double));
  }

  /* ================================================================== */
  /* Iterative Sum-Product Algorithm (flooding schedule)                */
  /* ================================================================== */
  for (iter = 0; iter < max_iter; iter++) {

    /* ------------------------ Check node update ------------------- */
    for (i = 0; i < M; i++) {
      for (k = 0; k < deg_c[i]; k++) {

        int j_idx = check_node[i][k];

        double prod_sign = 1.0;
        double sum_spf_val = 0.0;

        /* For each neighbor variable node (except j_idx) */
        for (j = 0; j < deg_c[i]; j++) {
          int var = check_node[i][j];
          if (j != k) {
            double x = LLR[var] + u[i][var];
            prod_sign *= sign_val(x);
            sum_spf_val += spf(fabs(x));
          }
        }

        v[i][j_idx] = prod_sign * spf(sum_spf_val);
      }
    }

    /* ------------------------ Variable node update ---------------- */
    for (j = 0; j < N; j++) {
      for (k = 0; k < deg_v[j]; k++) {

        int i_idx = variable_node[j][k];
        double sum_v = 0.0;

        /* Sum messages from all checks except i_idx */
        for (i = 0; i < deg_v[j]; i++) {
          int cnode = variable_node[j][i];
          if (i != k) {
            sum_v += v[cnode][j];
          }
        }
        u[i_idx][j] = sum_v;
      }
    }

    /* ------------------------ Tentative decision ------------------ */
    for (j = 0; j < N; j++) {
      double sum = LLR[j];
      for (i = 0; i < deg_v[j]; i++) {
        sum += v[variable_node[j][i]][j];
      }
      ecc[j] = (sum >= 0.0) ? 1 : 0;
    }

    /* ------------------------ Parity check H·ecc^T ---------------- */
    int parity_ok = 1;
    for (i = 0; i < M; i++) {
      int parity = 0;
      for (k = 0; k < deg_c[i]; k++) {
        parity ^= ecc[check_node[i][k]];
      }
      if (parity != 0) {
        parity_ok = 0;
        break;
      }
    }

    /* Early stopping if all parity checks satisfied */
    if (parity_ok) {
      break;
    }
  }

  /* ------------------------------------------------------------------ */
  /* Extract information bits (systematic part)                          */
  /* Assumes: codeword layout = [parity bits (N-K) | info bits (K)]      */
  /* ------------------------------------------------------------------ */
  for (i = 0; i < K; i++) {
    inf[i] = ecc[i + (N - K)];
  }

  /* ------------------------------------------------------------------ */
  /* Memory cleanup                                                     */
  /* ------------------------------------------------------------------ */
  for (i = 0; i < M; i++) {
    free(check_node[i]);
  }
  free(check_node);

  for (j = 0; j < N; j++) {
    free(variable_node[j]);
  }
  free(variable_node);

  free(deg_c);
  free(deg_v);

  for (i = 0; i < M; i++) {
    free(u[i]);
    free(v[i]);
  }
  free(u);
  free(v);
}

/* ========================================================================== */
/* Bit-wise LLR Computation from Per-Symbol Likelihoods                       */
/* ========================================================================== */
/**
 * @brief Compute bit-wise LLRs from symbol-wise likelihoods p(y | x_k).
 *
 * This helper assumes:
 *   - Modulation alphabet size E (e.g., 2, 4, 8, 16, ...)
 *   - E is a power of 2, so log2(E) is an integer = number of bits per symbol
 *   - For each symbol index k ∈ [0, E−1], its bit-label is derived from the
 *     binary representation of k:
 *
 *         symbol_bits[k][b] = (k >> b) & 1
 *
 * Bit-wise LLR for bit position b at symbol index i:
 *
 *     LLR[b + i * logE] = log( P(bit_b = 1 | y_i) / P(bit_b = 0 | y_i) )
 *                       ≈ log( Σ_{k:bit_b=1} p(y_i | x_k)
 *                              / Σ_{k:bit_b=0} p(y_i | x_k) )
 *
 * @param pyx  2D array pyx[E][N]: likelihood p(y_i | x_k) per symbol and time
 * @param E    Modulation order (number of symbols, must be power of 2)
 * @param N    Number of symbol positions (time index)
 * @param LLR  Output bit-wise LLRs, length = N * log2(E)
 */
void compute_llr_from_pyx(double **pyx, int E, int N, double *LLR) {
  int i, k, b;
  int logE =
      (int)log2((double)E); /* bits per symbol, assumes E is power of 2 */

  /* ------------------------------------------------------ */
  /* Precompute bit labels for each symbol index k          */
  /* ------------------------------------------------------ */
  int **symbol_bits = (int **)malloc(E * sizeof(int *));
  for (k = 0; k < E; k++) {
    symbol_bits[k] = (int *)malloc(logE * sizeof(int));
    for (b = 0; b < logE; b++) {
      symbol_bits[k][b] = (k >> b) & 1;
    }
  }

  /* ------------------------------------------------------ */
  /* Compute bit-wise LLR for each bit position and symbol  */
  /* ------------------------------------------------------ */
  for (i = 0; i < N; i++) {
    for (b = 0; b < logE; b++) {

      double p1 = 0.0;
      double p0 = 0.0;

      for (k = 0; k < E; k++) {
        if (symbol_bits[k][b])
          p1 += pyx[k][i];
        else
          p0 += pyx[k][i];
      }

      /* Numerical safety: avoid log(0) */
      if (p1 <= 0.0)
        p1 = 1e-300;
      if (p0 <= 0.0)
        p0 = 1e-300;

      LLR[b + i * logE] = log(p1 / p0);
    }
  }

  /* cleanup */
  for (k = 0; k < E; k++) {
    free(symbol_bits[k]);
  }
  free(symbol_bits);
}
