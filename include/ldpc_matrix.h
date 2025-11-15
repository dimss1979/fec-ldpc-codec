#ifndef LDPC_MATRIX_H
#define LDPC_MATRIX_H

#include <stdio.h>
#include <stdlib.h>

// Gallager法により (N_LDPC*w_c/w_r) × (N_LDPC) の検査行列を生成
void generate_Hmatrix(int **H_LDPC, int N_LDPC, int w_c, int w_r);

// 検査行列から (N_LDPC - N_LDPC*w_c/w_r) × (N_LDPC) の生成行列を生成
void generate_Gmatrix(int **H_LDPC, int **G_LDPC, int N_LDPC, int w_c, int w_r);

// 検査行列の4ループ数を計算
int count_floop(int **H_LDPC, int N_LDPC, int w_c, int w_r);

#endif // LDPC_MATRIX_H
