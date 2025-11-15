#include "ldpc_matrix.h"

// Gallager法により(N_LDPC*w_c/w_r)×(N_LDPC)の検査行列を生成
void generate_Hmatrix(int **H_LDPC, int N_LDPC, int w_c, int w_r) {
  // ローカル変数の定義
  int i, j, k, l, m;          // カウント用
  int buf;                    // 一時保存用
  int M = N_LDPC * w_c / w_r; // 検査行列の行数
  int N = N_LDPC;             // 検査行列の列数
  int B_M = M / w_c; // 検査行列のサブブロックの行数（小数点以下は切り捨て）

  // malloc関数（宣言）
  int *perm; // 順列
  perm = (int *)malloc(N * sizeof(int));

  // クリーニング
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      H_LDPC[i][j] = 0;
    }
  }

  // Gallager準備
  for (i = 0; i < B_M; i++) {
    for (j = i * w_r; j < w_r * (i + 1); j++) {
      H_LDPC[i][j] = 1;
    }
  }

  // Gallager法
  for (i = 1; i < w_c; i++) {
    // 整数を順番に格納
    for (l = 0; l < N; l++) {
      perm[l] = l;
    }
    // 列のランダム入れ替え
    for (l = 0; l < N; l++) {
      m = rand() % N;
      buf = perm[l];
      perm[l] = perm[m];
      perm[m] = buf;
    }
    for (j = B_M * i; j < B_M * (i + 1); j++) {
      for (k = 0; k < N; k++) {
        H_LDPC[j][k] = H_LDPC[j - B_M * i][perm[k]];
      }
    }
  }

  // malloc関数（解放）
  free(perm);
}

//(N_LDPC*w_c/w_r)×(N_LDPC)の検査行列から(N_LDPC-N_LDPC*w_c/w_r)×(N_LDPC)の生成行列を生成
void generate_Gmatrix(int **H_LDPC, int **G_LDPC, int N_LDPC, int w_c,
                      int w_r) {
  // ローカル変数の定義
  int i, j, k, l;             // カウント用
  int M = N_LDPC * w_c / w_r; // 検査行列の行数
  int N = N_LDPC;             // 検査行列の列数

  // malloc関数（宣言）
  int **X; // 変形用行列X
  X = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++)
    X[i] = (int *)malloc((M + N) * sizeof(int));

  int *X_rc; // 行コピー用
  X_rc = (int *)malloc((M + N) * sizeof(int));

  int *X_cc; // 列コピー用
  X_cc = (int *)malloc(N * sizeof(int));

  int *H_cc; // 列コピー用
  H_cc = (int *)malloc(M * sizeof(int));

  // 検査行列H_LDPCの転置行列H^Tを変形行列Xの左側に格納
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      X[i][j] = H_LDPC[j][i];
    }
  }

  // 変形行列Xの右側にN×N単位行列を格納
  for (i = 0; i < N; i++) {
    for (j = M; j < (M + N); j++) {
      if (i == j - M)
        X[i][j] = 1;
      else
        X[i][j] = 0;
    }
  }

  // 変形行列Xの左側の変形（ここの列基本変形はHに影響しない）
  for (j = 0; j < M; j++) {
    if (X[j][j] == 1) {
      for (i = 0; i < N; i++) {
        if (i != j && X[i][j] == 1) {
          for (k = 0; k < (M + N); k++) {
            X[i][k] = X[i][k] ^ X[j][k]; // 排他的論理和
          }
        }
      }
    } else {
      for (i = j + 1; i < N; i++) {
        if (X[i][j] == 1) {
          for (k = 0; k < (M + N); k++) {
            X_rc[k] = X[i][k]; // i行のコピー
            X[i][k] = X[j][k]; // 行交換j→i
            X[j][k] = X_rc[k]; // 行交換i→j
          }
          break;
        } else if (i == N - 1) {
          for (k = (M + N) - 1; k > j; k--) {
            if (X[j][k] == 1) {
              for (l = 0; l < N; l++) {
                X_cc[l] = X[l][k]; // 列のコピー
                X[l][k] = X[l][j]; // 列交換j→k
                X[l][j] = X_cc[l]; // 列交換k→j
              }
              break;
            }
          }
        }
      }
      for (i = 0; i < N; i++) {
        if (i != j && X[i][j] == 1) {
          for (k = 0; k < (M + N); k++) {
            X[i][k] = X[i][k] ^ X[j][k]; // 排他的論理和
          }
        }
      }
    }
  }

  // 変形行列Xの右側の変形（ここの列基本変形はHに影響する）
  for (j = 2 * M; j < (M + N); j++) {
    if (X[j - M][j] == 1) {
      for (i = 0; i < N; i++) {
        if (i != (j - M) && X[i][j] == 1) {
          for (k = 0; k < (M + N); k++) {
            X[i][k] = X[i][k] ^ X[j - M][k]; // 排他的論理和
          }
        }
      }
    } else {
      for (i = j - M; i < N; i++) {
        if (X[i][j] == 1) {
          for (k = 0; k < (M + N); k++) {
            X_rc[k] = X[i][k];     // i行のコピー
            X[i][k] = X[j - M][k]; // 行交換(j-M)→i
            X[j - M][k] = X_rc[k]; // 行交換i→(j-M)
          }
          break;
        } else if (i == N - 1) {
          for (k = (M + N) - 1; k > M - 1; k--) {
            if (X[j - M][k] == 1) {
              for (l = 0; l < N; l++) {
                X_cc[l] = X[l][k]; // 列のコピー
                X[l][k] = X[l][j]; // 列交換j→k
                X[l][j] = X_cc[l]; // 列交換k→j
              }
              for (l = 0; l < M; l++) {
                H_cc[l] = H_LDPC[l][k - M];          // 列のコピー
                H_LDPC[l][k - M] = H_LDPC[l][j - M]; // 列交換(j-M)→k
                H_LDPC[l][j - M] = H_cc[l];          // 列交換k→(j-M)
              } // 検査行列H_LDPCもXの右側と同様に列交換を行う
              break;
            }
          }
        }
      }
      for (i = 0; i < N; i++) {
        if (i != (j - M) && X[i][j] == 1) {
          for (k = 0; k < (M + N); k++) {
            X[i][k] = X[i][k] ^ X[j - M][k]; // 排他的論理和
          }
        }
      }
    }
  }

  // 生成行列G_LDPCの取り出し
  for (i = M; i < N; i++) {
    for (j = M; j < (M + N); j++) {
      G_LDPC[i - M][j - M] = X[i][j];
    }
  }

  // malloc関数（解放）
  for (i = 0; i < N; i++)
    free(X[i]);
  free(X);
  free(X_rc);
  free(X_cc);
  free(H_cc);
}

// nの階乗
static int factorial(int n) {
  // ローカル変数の定義
  int i;      // カウント用
  int fn = 1; // 階乗

  for (i = 0; i <= n; i++) {
    if (i == 0)
      fn = 1;
    else
      fn *= i;
  }

  return fn;
}

// 検査行列の4ループ数を計算する
int count_floop(int **H_LDPC, int N_LDPC, int w_c, int w_r) {
  // ローカル変数の定義
  int i, j, k, l;             // カウント用
  int M = N_LDPC * w_c / w_r; // 検査行列H_LDPCの行数
  int N = N_LDPC;             // 検査行列H_LDPCの列
  int mt;                     // 一致数
  int fl;                     // 4ループ数

  // malloc関数（宣言）
  int **vn; // 変数ノード
  vn = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++)
    vn[i] = (int *)malloc(w_c * sizeof(int));

  // 変数ノードの作成
  for (j = 0; j < N; j++) {
    k = 0;
    for (i = 0; i < M; i++) {
      if (H_LDPC[i][j] == 1) {
        vn[j][k] = i;
        k++;
        if (k == w_c)
          break; // forから出る
      }
    }
  }

  // 4ループ数の計算
  fl = 0;
  for (i = 0; i < N - 1; i++) {
    for (j = i + 1; j < N; j++) {
      mt = 0;
      for (k = 0; k < w_c; k++) {
        for (l = 0; l < w_c; l++) {
          if (vn[i][k] == vn[j][l])
            mt++;
        }
      }
      if (mt <= 1)
        fl += 0;
      else if (mt == 2)
        fl += 1;
      else
        fl +=
            factorial(mt) / factorial(mt - 2) / 2; // mt個から2個選ぶ組み合わせ
    }
  }

  // malloc関数（解放）
  for (i = 0; i < N; i++)
    free(vn[i]);
  free(vn);

  return fl;
}
