#ifndef LDPC_ENCODER_H
#define LDPC_ENCODER_H

// --------------------------------------------
// Function prototypes
// --------------------------------------------

// LDPC 検査行列 H と 生成行列 G の読み込み
void init_LDPC(int **H, int **G, int N_LDPC, int w_c, int w_r);

// LDPC 符号化 ecc = inf × G
void encode_LDPC(int *inf, int *ecc, int **G, int N_LDPC, int K_LDCP);

#endif
