/**
 * @file gene_hg.c
 * @brief LDPC H/G matrix generator (Gallager construction + G from H).
 *
 * This tool:
 *   1. Generates an LDPC parity-check matrix H via Gallager's regular
 * construction
 *   2. Constructs a systematic generator matrix G from H (GF(2) Gaussian
 * elimination)
 *   3. Counts 4-cycles in H (short cycles in the Tanner graph)
 *   4. Searches for the H/G pair with the smallest number of 4-cycles
 *   5. Periodically saves the best matrices and statistics into files
 *
 * Notes:
 *   - The search is performed by repeated random Gallager constructions.
 *   - For large N, the exhaustive search with a huge loop_count_max is
 *     computationally very expensive. Adjust loop_count_max as needed
 *     for practical use.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> /* mkdir() for POSIX */
#include <sys/types.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h> /* _mkdir() on Windows */
#endif

#include "ldpc_matrix.h"

/* ------------------------------------------------------------------------- */
/* Portable mkdir wrapper                                                    */
/*   - On Windows, uses _mkdir()                                            */
/*   - On POSIX systems, uses mkdir(path, mode)                             */
/* ------------------------------------------------------------------------- */
static void make_dir(const char *d) {
#ifdef _WIN32
  _mkdir(d);
#else
  mkdir(d, 0755);
#endif
}

/* ========================================================================== */
/* MAIN                                                                       */
/* ========================================================================== */
int main(void) {
  srand((unsigned int)time(NULL));

  printf("==============================================\n");
  printf("       LDPC Matrix Generator (Gallager)       \n");
  printf("==============================================\n\n");

  /* ------------------------------------------------------------------ */
  /* User input: (N, wc, wr)                                            */
  /* ------------------------------------------------------------------ */
  int N, wc, wr;

  printf("Codeword length N: ");
  scanf("%d", &N);

  printf("Column weight wc (small: 2 or 3): ");
  scanf("%d", &wc);

  printf("Row weight wr (larger than wc): ");
  scanf("%d", &wr);

  int M = (N * wc) / wr; /* number of parity-check equations */
  int K = N - M;         /* number of information bits       */
  double R = (double)K / N;

  printf("\nRate R = %.5f (K = %d, M = %d)\n\n", R, K, M);

  /* ------------------------------------------------------------------ */
  /* Prepare output directory                                           */
  /*   matrices/Nxxx_wcX_wrY/                                          */
  /* ------------------------------------------------------------------ */
  char dirpath[128];
  sprintf(dirpath, "matrices/N%d_wc%d_wr%d", N, wc, wr);

  make_dir("matrices");
  make_dir(dirpath);

  /* Output file paths */
  char path_H[256], path_G[256], path_info[256];
  sprintf(path_H, "%s/H.csv", dirpath);
  sprintf(path_G, "%s/G.csv", dirpath);
  sprintf(path_info, "%s/info.txt", dirpath);

  /* ------------------------------------------------------------------ */
  /* Allocate matrices H, G and their "best" copies                     */
  /* ------------------------------------------------------------------ */
  int **H = (int **)malloc(M * sizeof(int *));
  int **H_best = (int **)malloc(M * sizeof(int *));
  for (int i = 0; i < M; i++) {
    H[i] = (int *)malloc(N * sizeof(int));
    H_best[i] = (int *)malloc(N * sizeof(int));
  }

  int **G = (int **)malloc(K * sizeof(int *));
  int **G_best = (int **)malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    G[i] = (int *)malloc(N * sizeof(int));
    G_best[i] = (int *)malloc(N * sizeof(int));
  }

  /* ------------------------------------------------------------------ */
  /* Search H/G matrices with minimum number of 4-cycles                */
  /* ------------------------------------------------------------------ */
  const int loop_count_max = 10000000;
  const double print_interval_sec = 1.0; /* periodic save interval   */

  int loop;
  int best_floop = -1;     /* best (minimum) number of 4-cycles */
  long long floop_sum = 0; /* average tracking */

  clock_t t_start = clock();
  clock_t t_last_print = clock();

  printf("Searching for best H/G matrices (min 4-cycles)...\n");

  for (loop = 1; loop <= loop_count_max; loop++) {

    /* -------------------------------------------------------------- */
    /* 1) Generate new H and G                                        */
    /* -------------------------------------------------------------- */
    generate_Hmatrix(H, N, wc, wr);
    generate_Gmatrix(H, G, N, wc, wr);

    /* -------------------------------------------------------------- */
    /* 2) Count 4-cycles in H                                         */
    /* -------------------------------------------------------------- */
    int floop = count_floop(H, N, wc, wr);
    floop_sum += floop;

    /* -------------------------------------------------------------- */
    /* 3) Update best H/G if this one has fewer 4-cycles              */
    /* -------------------------------------------------------------- */
    if (best_floop == -1 || floop < best_floop) {
      best_floop = floop;

      /* copy best H */
      for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
          H_best[i][j] = H[i][j];

      /* copy best G */
      for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
          G_best[i][j] = G[i][j];
    }

    /* -------------------------------------------------------------- */
    /* 4) Periodically save current best matrices and statistics       */
    /* -------------------------------------------------------------- */
    clock_t t_now = clock();
    double elapsed = (double)(t_now - t_last_print) / CLOCKS_PER_SEC;

    if (loop == 1 || elapsed > print_interval_sec) {

      t_last_print = t_now;

      /* Save best H as CSV (no separators, 0/1 matrix) */
      FILE *fp = fopen(path_H, "w");
      if (fp) {
        for (int i = 0; i < M; i++) {
          for (int j = 0; j < N; j++)
            fprintf(fp, "%d", H_best[i][j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }

      /* Save best G as CSV (no separators, 0/1 matrix) */
      fp = fopen(path_G, "w");
      if (fp) {
        for (int i = 0; i < K; i++) {
          for (int j = 0; j < N; j++)
            fprintf(fp, "%d", G_best[i][j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }

      /* Save status information */
      fp = fopen(path_info, "w");
      if (fp) {
        fprintf(fp, "LDPC Matrix Generation Status\n");
        fprintf(fp, "Code rate R = %.5f\n", R);
        fprintf(fp, "N = %d\n", N);
        fprintf(fp, "wc = %d\n", wc);
        fprintf(fp, "wr = %d\n", wr);
        fprintf(fp, "Loop count = %d\n", loop);
        fprintf(fp, "Best 4-cycles = %d\n", best_floop);
        fprintf(fp, "Average 4-cycles = %.3f\n", (double)floop_sum / loop);
        fclose(fp);
      }

      printf("[Loop %d] Best 4-cycles = %d, Avg = %.3f\n", loop, best_floop,
             (double)floop_sum / loop);
    }
  }

  /* ------------------------------------------------------------------ */
  /* Cleanup                                                            */
  /* ------------------------------------------------------------------ */
  for (int i = 0; i < M; i++) {
    free(H[i]);
    free(H_best[i]);
  }
  free(H);
  free(H_best);

  for (int i = 0; i < K; i++) {
    free(G[i]);
    free(G_best[i]);
  }
  free(G);
  free(G_best);

  printf("\nGeneration completed.\n");
  printf("Files saved under directory: %s\n", dirpath);

  (void)t_start; /* currently unused; kept for potential profiling */

  return 0;
}
