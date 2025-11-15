#include "ldpc_matrix.h"
#include <direct.h> // Windowsでディレクトリ作成用 (_mkdir)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
  srand((unsigned int)time(NULL));

  printf("=== LDPC行列生成プログラム ===\n");
  printf("Gallager法による検査行列・生成行列を生成し、CSVファイルに出力します。"
         "\n\n");

  int N_LDPC, w_c, w_r;
  int T_optitv = 1;

  printf("符号語長: ");
  scanf("%d", &N_LDPC);
  printf("検査行列H_LDPCの列重み（小さめ, 2, 3）: ");
  scanf("%d", &w_c);
  printf("検査行列H_LDPCの行重み（大きめ）: ");
  scanf("%d", &w_r);

  int M = N_LDPC * w_c / w_r;
  int N = N_LDPC;
  double R = 1 - (double)w_c / (double)w_r;

  // 出力フォルダ名を決定
  char dirname[128];
  sprintf(dirname, "matrices\\N%d_wc%d_wr%d", N_LDPC, w_c, w_r);

  // フォルダ作成（存在しない場合のみ）
  _mkdir("matrices"); // 親フォルダ
  _mkdir(dirname);    // サブフォルダ

  // 出力ファイル名
  char csvfilename_H[256], csvfilename_G[256], csvfilename_ch[256];
  sprintf(csvfilename_H, "%s\\H.csv", dirname);
  sprintf(csvfilename_G, "%s\\G.csv", dirname);
  sprintf(csvfilename_ch, "%s\\check.csv", dirname);

  printf("%s を実行\n", csvfilename_ch);

  // malloc関数（宣言）
  int **H_LDPC = (int **)malloc(M * sizeof(int *));
  int **H_min = (int **)malloc(M * sizeof(int *));
  for (int i = 0; i < M; i++) {
    H_LDPC[i] = (int *)malloc(N * sizeof(int));
    H_min[i] = (int *)malloc(N * sizeof(int));
  }

  int **G_LDPC = (int **)malloc((N - M) * sizeof(int *));
  int **G_min = (int **)malloc((N - M) * sizeof(int *));
  for (int i = 0; i < (N - M); i++) {
    G_LDPC[i] = (int *)malloc(N * sizeof(int));
    G_min[i] = (int *)malloc(N * sizeof(int));
  }

  // ループ関係の変数
  int loop = 1;
  int loop_max = 2000000000; // 2*10^9
  int fl1 = 0, fl2 = 0;
  double fl_avg = 0;
  double prog = 0, totm = 0;
  clock_t start, end;
  FILE *fp;

  // ループ開始準備
  start = clock();
  for (loop = 1; loop < loop_max; loop++) {
    generate_Hmatrix(H_LDPC, N_LDPC, w_c, w_r);
    generate_Gmatrix(H_LDPC, G_LDPC, N_LDPC, w_c, w_r);

    // 4ループ数の計算と評価
    fl1 = count_floop(H_LDPC, N_LDPC, w_c, w_r);
    fl_avg += fl1;
    if (loop == 1 || fl1 < fl2) {
      fl2 = fl1;
      for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
          H_min[i][j] = H_LDPC[i][j];
      for (int i = 0; i < N - M; i++)
        for (int j = 0; j < N; j++)
          G_min[i][j] = G_LDPC[i][j];
    }

    end = clock();
    prog = (double)(end - start) / CLOCKS_PER_SEC;
    if (loop == 1 || prog > T_optitv) {
      totm += prog;

      fp = fopen(csvfilename_H, "w");
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++)
          fprintf(fp, "%d", H_min[i][j]);
        fprintf(fp, "\n");
      }
      fclose(fp);

      fp = fopen(csvfilename_G, "w");
      for (int i = 0; i < N - M; i++) {
        for (int j = 0; j < N; j++)
          fprintf(fp, "%d", G_min[i][j]);
        fprintf(fp, "\n");
      }
      fclose(fp);

      fp = fopen(csvfilename_ch, "w");
      fprintf(fp, "検査行列と生成行列の算出状況\n");
      fprintf(fp, "符号化率 %.5f\n", R);
      fprintf(fp, "符号語ビット長 %d\n", N_LDPC);
      fprintf(fp, "列重み %d\n", w_c);
      fprintf(fp, "行重み %d\n", w_r);
      fprintf(fp, "生成回数 %d\n", loop);
      fprintf(fp, "生成時間 %.0f\n", totm);
      fprintf(fp, "現時点の4ループ数 %d\n", fl2);
      fprintf(fp, "4ループ数の平均値 %.2f\n", (double)fl_avg / loop);
      fclose(fp);

      start = clock();
    }
  }

  // malloc関数（解放）
  for (int i = 0; i < M; i++) {
    free(H_LDPC[i]);
    free(H_min[i]);
  }
  free(H_LDPC);
  free(H_min);

  for (int i = 0; i < N - M; i++) {
    free(G_LDPC[i]);
    free(G_min[i]);
  }
  free(G_LDPC);
  free(G_min);

  printf("\n処理が完了しました。生成されたCSVファイルを確認してください。\n");
  return 0;
}
