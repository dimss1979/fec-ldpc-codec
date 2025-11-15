#include "ldpc_encoder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void init_LDPC(int **H, int **G, int N_LDPC, int w_c, int w_r) {
  int i, j;
  char path_H[256];
  char path_G[256];
  char line[4096];
  FILE *fp;

  // ----------------------------------------------------
  // パス生成： matrices/N1024_wc3_wr6/H.csv
  // ----------------------------------------------------
  sprintf(path_H, "matrices/N%d_wc%d_wr%d/H.csv", N_LDPC, w_c, w_r);
  sprintf(path_G, "matrices/N%d_wc%d_wr%d/G.csv", N_LDPC, w_c, w_r);

  // ----------------------------------------------------
  // H の読み込み
  // ----------------------------------------------------
  fp = fopen(path_H, "r");
  if (!fp) {
    printf("エラー: %s が見つかりません\n", path_H);
    exit(1);
  }

  i = 0;
  while (fgets(line, sizeof(line), fp)) {
    for (j = 0; j < N_LDPC; j++) {
      H[i][j] = (line[j] == '1') ? 1 : 0;
    }
    i++;
  }
  fclose(fp);

  // ----------------------------------------------------
  // G の読み込み
  // ----------------------------------------------------
  fp = fopen(path_G, "r");
  if (!fp) {
    printf("エラー: %s が見つかりません\n", path_G);
    exit(1);
  }

  i = 0;
  while (fgets(line, sizeof(line), fp)) {
    for (j = 0; j < N_LDPC; j++) {
      G[i][j] = (line[j] == '1') ? 1 : 0;
    }
    i++;
  }
  fclose(fp);
}

// ============================================================
// LDPC エンコード ecc = inf × G (mod 2)
// ------------------------------------------------------------
// inf  : 情報ビット (長さ K_LDPC)
// ecc  : 符号語 (長さ N_LDPC)
// G    : K×N 生成行列
// ============================================================
void encode_LDPC(int *inf, int *ecc, int **G, int N_LDPC, int K_LDPC) {
  int i, j;

  for (i = 0; i < N_LDPC; i++) {
    int val = 0;
    for (j = 0; j < K_LDPC; j++) {
      val ^= (inf[j] & G[j][i]); // XOR accumulation
    }
    ecc[i] = val;
  }
}
