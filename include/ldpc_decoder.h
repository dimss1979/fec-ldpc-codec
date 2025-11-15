#ifndef LDPC_DECODER_H
#define LDPC_DECODER_H

// ------------------------------------
// Sum-Product LDPC Decoder (SPA)
// ------------------------------------
// Inputs:
//   LLR[N]    : Received LLR from channel
//   H         : Parity check matrix (M × N)
//   M         : Number of parity equations = N - K
//   N         : Codeword length
//   K         : Number of information bits
//   max_iter  : Maximum SPA iterations
//
// Outputs:
//   ecc[N]    : Decoded codeword
//   inf[K]    : Decoded information bits
//-------------------------------------

void ldpc_decode_spa(double *LLR, int *ecc, int *inf, int **H, int M, int N,
                     int K, int max_iter);

// ------------------------------------
// Compute LLR from likelihood pyx
// pyx[E][N] : symbol likelihood
// E         : modulation order (e.g., E=4 for QPSK)
// N         : number of bits
//------------------------------------
void compute_llr_from_pyx(double **pyx, int E, int N, double *LLR);

#endif
