/* Link-time DGEMM call-shape profiler for the wp=8 benchmark build.
 *
 * Link an executable with this object and -Wl,--wrap=dgemm_.  The wrapper
 * records scalar arguments without changing the bundled implementation, then
 * prints one aggregate record per distinct call shape at normal process exit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { max_shapes = 512 };

struct dgemm_shape {
    char transa;
    char transb;
    int m;
    int n;
    int k;
    int lda;
    int ldb;
    int ldc;
    double alpha;
    double beta;
    unsigned long long calls;
    unsigned long long multiply_adds;
};

static struct dgemm_shape shapes[max_shapes];
static int shape_count;
static int report_registered;

extern void __real_dgemm_(char *, char *, int *, int *, int *, double *,
                          double *, int *, double *, int *, double *,
                          double *, int *, size_t, size_t);

static void report_shapes(void)
{
    int index;

    for (index = 0; index < shape_count; ++index) {
        const struct dgemm_shape *shape = &shapes[index];
        fprintf(stderr,
                "DGEMM_SHAPE,%c,%c,%d,%d,%d,%.17g,%.17g,%d,%d,%d,%llu,%llu\n",
                shape->transa, shape->transb, shape->m, shape->n, shape->k,
                shape->alpha, shape->beta, shape->lda, shape->ldb, shape->ldc,
                shape->calls, shape->multiply_adds);
    }
}

static void record_shape(char transa, char transb, int m, int n, int k,
                         double alpha, double beta, int lda, int ldb, int ldc)
{
    int index;

    for (index = 0; index < shape_count; ++index) {
        struct dgemm_shape *shape = &shapes[index];
        if (shape->transa == transa && shape->transb == transb &&
            shape->m == m && shape->n == n && shape->k == k &&
            shape->lda == lda && shape->ldb == ldb && shape->ldc == ldc &&
            shape->alpha == alpha && shape->beta == beta) {
            ++shape->calls;
            shape->multiply_adds +=
                (unsigned long long)m * (unsigned long long)n *
                (unsigned long long)k;
            return;
        }
    }

    if (shape_count >= max_shapes) {
        fputs("DGEMM profile exceeded max_shapes\n", stderr);
        exit(EXIT_FAILURE);
    }
    shapes[shape_count] = (struct dgemm_shape){
        transa, transb, m, n, k, lda, ldb, ldc, alpha, beta, 1,
        (unsigned long long)m * (unsigned long long)n * (unsigned long long)k
    };
    ++shape_count;
}

void __wrap_dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc,
                   size_t transa_length, size_t transb_length)
{
    if (!report_registered) {
        if (atexit(report_shapes) != 0) {
            fputs("Cannot register DGEMM profile reporter\n", stderr);
            exit(EXIT_FAILURE);
        }
        report_registered = 1;
    }
    record_shape(*transa, *transb, *m, *n, *k, *alpha, *beta,
                 *lda, *ldb, *ldc);
    __real_dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb,
                  beta, c, ldc, transa_length, transb_length);
}
