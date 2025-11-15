#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ldpc_decoder.h"
#include "ldpc_encoder.h"

#define PI 3.141592

/* ============================================================
 * Utility: 2D matrix allocation (int)
 * ============================================================ */
static int **alloc_matrix_int(int rows, int cols) {
  int **m = (int **)malloc(rows * sizeof(int *));
  if (!m) {
    fprintf(stderr, "alloc_matrix_int: malloc failed (rows)\n");
    exit(1);
  }
  for (int i = 0; i < rows; i++) {
    m[i] = (int *)malloc(cols * sizeof(int));
    if (!m[i]) {
      fprintf(stderr, "alloc_matrix_int: malloc failed (row %d)\n", i);
      exit(1);
    }
  }
  return m;
}

static void free_matrix_int(int **m, int rows) {
  for (int i = 0; i < rows; i++) {
    free(m[i]);
  }
  free(m);
}

/* ============================================================
 * Utility: Gaussian random (Box-Muller)
 * ============================================================ */
static double rand_uniform(void) { return (rand() + 1.0) / (RAND_MAX + 2.0); }

static double randn(void) {
  double u1 = rand_uniform();
  double u2 = rand_uniform();
  return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ============================================================
 * Main: LDPC BER over AWGN (BPSK)
 * ============================================================ */
int main(void) {
  /* --------------------------------------------------------
   * 1. LDPC parameters (適宜変更可能)
   * -------------------------------------------------------- */
  int N_LDPC = 1024;
  int w_c = 3;
  int w_r = 6;

  /* M: parity bits, K: information bits */
  int M_LDPC = (N_LDPC * w_c) / w_r;
  int K_LDPC = N_LDPC - M_LDPC;

  printf("LDPC Parameters:\n");
  printf("  N = %d\n", N_LDPC);
  printf("  K = %d\n", K_LDPC);
  printf("  M = %d\n", M_LDPC);
  printf("  wc = %d, wr = %d\n\n", w_c, w_r);

  /* --------------------------------------------------------
   * 2. Matrices H (M×N), G (K×N) を読み込み
   * -------------------------------------------------------- */
  int **H = alloc_matrix_int(M_LDPC, N_LDPC);
  int **G = alloc_matrix_int(K_LDPC, N_LDPC);

  printf("Loading H and G... ");
  init_LDPC(H, G, N_LDPC, w_c, w_r);
  printf("Done.\n\n");

  /* --------------------------------------------------------
   * 3. Simulation settings
   * -------------------------------------------------------- */
  const double EbN0_dB_list[] = {0, 1, 2, 3, 4, 5};
  const int NUM_SNR = (int)(sizeof(EbN0_dB_list) / sizeof(EbN0_dB_list[0]));

  int frames_per_snr = 200; /* フレーム数（必要に応じて増やす） */
  int max_iter = 30;        /* SPA 最大反復回数 */

  /* 乱数シード */
  srand((unsigned int)time(NULL));

  /* 作業用バッファ */
  int *inf = (int *)malloc(K_LDPC * sizeof(int));         /* 送信情報ビット */
  int *code = (int *)malloc(N_LDPC * sizeof(int));        /* 送信符号語ビット */
  double *tx = (double *)malloc(N_LDPC * sizeof(double)); /* BPSK symbols */
  double *rx = (double *)malloc(N_LDPC * sizeof(double)); /* 受信値 */
  double *LLR = (double *)malloc(N_LDPC * sizeof(double)); /* 受信LLR */
  int *ecc_hat = (int *)malloc(N_LDPC * sizeof(int)); /* 復号後符号語ビット */
  int *inf_hat = (int *)malloc(K_LDPC * sizeof(int)); /* 復号後情報ビット */

  if (!inf || !code || !tx || !rx || !LLR || !ecc_hat || !inf_hat) {
    fprintf(stderr, "malloc failed.\n");
    return 1;
  }

  /* 結果CSVの出力ファイル */
  FILE *fp_csv = fopen("ldpc_ber_awgn.csv", "w");
  if (!fp_csv) {
    fprintf(stderr, "Cannot open ldpc_ber_awgn.csv for write.\n");
    return 1;
  }
  fprintf(fp_csv, "EbN0_dB,BER_info,BER_code\n");

  /* --------------------------------------------------------
   * 4. SNR loop
   * -------------------------------------------------------- */
  for (int idx_snr = 0; idx_snr < NUM_SNR; idx_snr++) {
    double EbN0_dB = EbN0_dB_list[idx_snr];
    double EbN0 = pow(10.0, EbN0_dB / 10.0);

    double R = (double)K_LDPC / (double)N_LDPC; /* 符号化率 */
    /* BPSK: Es=1 として、sigma^2 = N0/2,  Eb = Es/R → N0 = Eb/EbN0 → sigma^2 =
     * Eb/(2EbN0) = 1/(2 R EbN0) */
    double sigma2 = 1.0 / (2.0 * R * EbN0);
    double sigma = sqrt(sigma2);

    long long err_info = 0;
    long long err_code = 0;
    long long total_info_bits = (long long)frames_per_snr * K_LDPC;
    long long total_code_bits = (long long)frames_per_snr * N_LDPC;

    printf("SNR = %.2f dB (R=%.3f, sigma^2=%.4g)\n", EbN0_dB, R, sigma2);

    for (int frm = 0; frm < frames_per_snr; frm++) {

      /* ---------------------------------------------
       * (1) 情報ビット生成（ランダム）
       * --------------------------------------------- */
      for (int i = 0; i < K_LDPC; i++) {
        inf[i] = rand() & 1;
      }

      /* ---------------------------------------------
       * (2) LDPC 符号化
       * code[]: Nビットの符号語
       * --------------------------------------------- */
      encode_LDPC(inf, code, G, N_LDPC, K_LDPC);

      /* ---------------------------------------------
       * (3) BPSK 変調: 0→-1, 1→+1
       * --------------------------------------------- */
      for (int i = 0; i < N_LDPC; i++) {
        tx[i] = (code[i] == 1) ? +1.0 : -1.0;
      }

      /* ---------------------------------------------
       * (4) AWGN チャネル
       * --------------------------------------------- */
      for (int i = 0; i < N_LDPC; i++) {
        double n = sigma * randn();
        rx[i] = tx[i] + n;
      }

      /* ---------------------------------------------
       * (5) LLR 計算
       *    mapping: bit=1 → +1, bit=0 → -1
       *    LLR = ln p(y|1)/p(y|0) = 2y / sigma^2
       * --------------------------------------------- */
      for (int i = 0; i < N_LDPC; i++) {
        LLR[i] = 2.0 * rx[i] / sigma2;
      }

      /* ---------------------------------------------
       * (6) Sum-Product 復号
       *    ecc_hat: Nビット符号語
       *    inf_hat: Kビット情報
       * --------------------------------------------- */
      ldpc_decode_spa(LLR, ecc_hat, inf_hat, H, M_LDPC, N_LDPC, K_LDPC,
                      max_iter);

      /* ---------------------------------------------
       * (7) 誤りビット数をカウント
       * --------------------------------------------- */
      for (int i = 0; i < K_LDPC; i++) {
        if (inf[i] != inf_hat[i])
          err_info++;
      }
      for (int i = 0; i < N_LDPC; i++) {
        if (code[i] != ecc_hat[i])
          err_code++;
      }
    }

    double BER_info = (double)err_info / (double)total_info_bits;
    double BER_code = (double)err_code / (double)total_code_bits;

    printf("  Info BER = %e  (errors=%lld / %lld)\n", BER_info, err_info,
           total_info_bits);
    printf("  Code BER = %e  (errors=%lld / %lld)\n\n", BER_code, err_code,
           total_code_bits);

    fprintf(fp_csv, "%.2f,%.10e,%.10e\n", EbN0_dB, BER_info, BER_code);
    fflush(fp_csv);
  }

  fclose(fp_csv);

  /* 後処理 */
  free(inf);
  free(code);
  free(tx);
  free(rx);
  free(LLR);
  free(ecc_hat);
  free(inf_hat);

  free_matrix_int(H, M_LDPC);
  free_matrix_int(G, K_LDPC);

  printf("LDPC BER simulation finished. Results saved to ldpc_ber_awgn.csv\n");

  return 0;
}
