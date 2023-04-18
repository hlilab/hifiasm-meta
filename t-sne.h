#ifndef T_SNE_H
#define T_SNE_H

typedef struct {
	int max_itr;       // maximum number of iterations
	double perplexity; // perplexity
	double theta;      // gradient accuracy
} ts_opt_t;

#ifdef __cplusplus
extern "C" {
#endif

double *ts_fit(int N, int n_in, double *x, int n_out, double theta, double perplexity, int seed);

double *sann_data_read_1d(const char *fn, int *n, int *n_col, char ***row_names, char ***col_names);
float **sann_data_read(const char *fn, int *n_, int *n_col_, char ***row_names, char ***col_names);
void sann_free_names(int n, char **s);
void sann_free_vectors(int n, float **x);

#ifdef __cplusplus
}
#endif

#endif
