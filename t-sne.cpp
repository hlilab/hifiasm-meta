#include <stdlib.h>
#include "tsne.h"
#include "t-sne.h"

double *ts_fit(int N, int n_in, double *x, int n_out, double theta, double perplexity, int seed)
{
	double *y;
	TSNE *ts = new TSNE();
	y = (double*)malloc(N * n_out * sizeof(double));
	ts->run(x, N, n_in, y, n_out, perplexity, theta, seed, false);
	delete(ts);
	return y;
}
