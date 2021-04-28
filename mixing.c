#include <stdio.h>
#define __USE_XOPEN
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>

#include <time.h>
#ifndef __unix__
#include <sys/time.h>
#endif

double MEPS = 1e-24;

#define NS_PER_SEC 1000000000
int64_t wall_clock_ns()
{
#ifdef __unix__
	struct timespec tspec;
	int r = clock_gettime(CLOCK_MONOTONIC, &tspec);
	assert(r==0);
	return tspec.tv_sec*NS_PER_SEC + tspec.tv_nsec;
#else
	struct timeval tv;
	int r = gettimeofday( &tv, NULL );
	assert(r==0);
	return tv.tv_sec*NS_PER_SEC + tv.tv_usec*1000;
#endif
}

#define Malloc(size) malloc_or_abort((size), __LINE__)
void* malloc_or_abort(size_t size, int line)
{
    void *ptr = malloc(size);
    if(ptr == NULL) {
        fprintf(stderr, "Error not enough memory at %d (size=%zu)\n", line, size);
        exit(1);
    }
    return ptr;
}

#define Calloc(count, size) calloc_or_abort(count, (size), __LINE__)
void* calloc_or_abort(size_t count, size_t size, int line)
{
    void *ptr = calloc(count, size);
    if(ptr == NULL) {
        fprintf(stderr, "Errror not enough memory at %d (count=%zu, size=%zu)\n", line, count, size);
        exit(1);
    }
    return ptr;
}

double wall_time_diff(int64_t ed, int64_t st)
{
	return (double)(ed-st)/(double)NS_PER_SEC;
}

struct Node {
    int index;
    double value;
};
typedef struct Node Node;

struct CUTProblem {
    int n;
    int nnz;
    Node **C;
};
typedef struct CUTProblem CUTProblem;

struct Model {
    int n, k;
    double **V;
};
typedef struct Model Model;

enum {MAXCUT=0, MAXSAT};
struct Parameter {
    int solver;
	int k;
    double eps;
    int max_iter;
    int is_unspeficied_wcnf;
    int verbose;
    int n_trial;
    FILE *fin;
    char *fin_name;
};
typedef struct Parameter Parameter;

struct SATProblem {
    int m; // number of clauses
    int n; // number of variables
    int nnz; // number of literals
    int **index; // index[j] = [1 2 -4]
    int *clause_len, *var_len;
    int *clause_weight;
};
typedef struct SATProblem SATProblem;

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;
    // init
    if(max_line_len == 0) {
        max_line_len = 1024;
	    line = (char *)Malloc(max_line_len*sizeof(char));
    }

	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}

int is_int(char *s)
{
    if(!s || *s == '\0')
        return 0;
    if(*s!='-' && !isdigit(*s)) 
        return 0;
    s++;
    for(; *s; s++) {
        if(!isdigit(*s))
            return 0;
    }
    return 1;
}

int is_valid_literal(char *s, int nvar)
{
    if(!is_int(s)) return 0;
    int var = abs(atoi(s));
    if(var != 0 && var <= nvar)
        return 1;
    return 0;
}

// DIMAC format
// Error detection: If cnf and all clause start with 1, then the format is wcnf.
void read_sat_problem(FILE *fin, SATProblem *prob, int is_unspeficied_wcnf)
{
    char *delim = " \t\r\n";
    int has_weight = 0;
    int nnz=0, n=0, m=0;

    int *var_len = NULL;
    int *clause_len = NULL;

    int clause_cnt = 1;
    int lineno = 1;
    for(; readline(fin); lineno++) {
        if(line[0] == '\0' || strchr("c\r\n", line[0]))
            continue;
        if(line[0] == 'p') {
            char *p = strtok(line+1, delim);
            if(!strcmp(p, "cnf")) {
                has_weight = 0;
            } else if(!strcmp(p, "wcnf")) {
                has_weight = 1;
            } else {
                fprintf(stderr, "(line %d) No format specified (cnf or wcnf).\n", lineno);
                exit(1);
            }
            // some public dataset is wcnf but marked as cnf
            if(is_unspeficied_wcnf) has_weight = 1;

            p = strtok(NULL, delim);
            if(p)
                n = atoi(p);
            p = strtok(NULL, delim);
            if(p)
                m = atoi(p);
            if(m == 0 || n == 0) {
                fprintf(stderr, "(line %d) Wrong format in parameter line\n", lineno);
                exit(1);
            }
            clause_len = Calloc(m+1, sizeof(*clause_len));
            var_len = Calloc(n+1, sizeof(*var_len));
            continue;
        }
        if(var_len == NULL || clause_len == NULL) {
            fprintf(stderr, "(line %d) Clause appears before parameters\n", lineno);
            exit(1);
        }

        // count nnz and check has_weight
        int has_error = 0;
        char *p = strtok(line, delim);
        if(!p){
            fprintf(stderr, "(line %d) Empty line in clause\n", lineno);
            exit(1);
        }
        if(!strncmp(p, "1 %", 3)) // ending symbol in some format
            break;
        if(has_weight){
            if(!is_int(p) || *p=='-' || !strcmp(p, "0")){
                fprintf(stderr, "(line %d) Only accept positve integer weight\n", lineno);
                exit(1);
            }
            p = strtok(NULL, delim);
        }
        
        int is_zero_terminated = 0;

        for(; p; ){
            if(!strcmp(p, "0")){
                is_zero_terminated = 1;
                break;
            }
            nnz++;
            if(!is_valid_literal(p, n)){
                has_error = 1;
                break;
            }else{
                int var = abs(atoi(p));
                var_len[var]++;
                clause_len[clause_cnt]++;
            }
            p = strtok(NULL, delim);
        }
        if(clause_len[clause_cnt] == 0){
            fprintf(stderr, "(line %d) Clause has no literal\n", lineno);
            exit(1);
        }
        if(!is_zero_terminated){
            fprintf(stderr, "(line %d) Clause need to be terminated with 0\n", lineno);
        }
        if(has_error) {
            fprintf(stderr, "(line %d) Wrong format in clause\n", lineno);
            exit(1);
        }
        clause_cnt++;
        if(is_unspeficied_wcnf && clause_cnt == m+1) {
            fprintf(stderr, "Read all the clauses for irregular WCNF format. Ignored error check since line %d.\n", lineno);
            break;
        }else if(clause_cnt > m+1){
            fprintf(stderr, "(line %d) More clause then speficied in parameters (%d)\n", lineno, m);
            exit(1);
        }
    }
    if(clause_cnt < m+1){ // include the 0-th clause
        fprintf(stderr, "(error) Fewer clauses (%d) then speficied in parameters (%d) \n", clause_cnt, m);
        exit(1);
    }

    int *pool = (int*) Calloc(nnz+n + m+1, sizeof(*pool));
    int cur = 0;

    int **index = (int**) Calloc(n+1, sizeof(*index));
    int *clause_weight = (int*) Calloc(m+1, sizeof(*clause_weight));

    int *var_pos = (int*) Calloc(n+1, sizeof(*var_pos));
    var_pos[0] = 0;
    var_len[0] = m;
    for(int j=1; j<=n; j++) {
        var_pos[j] = var_pos[j-1] + var_len[j-1]+1;
        index[j] = pool + var_pos[j];
    }
#if 0
    printf("var_pos : ");
    for(int j=0; j<=n; j++)
        printf("%d ", var_pos[j]);
    printf("\n");
#endif

    // initialize truth vector
    index[0] = pool;
    for(int i=1; i<=m; i++) {
        pool[cur++] = i;
    }
    pool[cur++] = 0;

    // initialize variables
    fseek(fin, 0, SEEK_SET);

    clause_cnt = 1;
    lineno = 1;
    for(; readline(fin); lineno++) {
        if(line[0] == '\0' || strchr("c\n", line[0]))
            continue;
        if(line[0] == 'p')
            continue;
        char *p = strtok(line, delim);
        if(has_weight){
            clause_weight[clause_cnt] = atoi(p);
            p = strtok(NULL, delim);
        }else{
            clause_weight[clause_cnt] = 1;
        }
        while(p){
            int literal = atoi(p);
            if(literal == 0)
                break;
            int sign = (literal>0)? 1 : -1;
            int var = abs(literal);
            //printf("lit %d at pos %d = %d\n", literal, var_pos[var], sign*clause_cnt);
            pool[var_pos[var]++] = sign * clause_cnt;

            p = strtok(NULL, delim);
        }
        clause_cnt++;
    }

    prob->m = m;
    prob->n = n;
    prob->nnz = nnz;
    prob->clause_len = clause_len;
    prob->var_len = var_len;
    prob->index = index;
    prob->clause_weight = clause_weight;
    free(var_pos);
}

void print_struct(SATProblem *prob)
{
    printf("n %d m %d nnz %d\n", prob->n, prob->m, prob->nnz);
    printf("clause_len: ");
    for(int i=1; i<=prob->m; i++)
        printf("%d ", prob->clause_len[i]);
    printf("\n var:\n");
    for(int j=1; j<=prob->n; j++){
        printf("len %d : ", prob->var_len[j]);
        int *p = prob->index[j];
        for(; *p; p++){
            printf("%d ", *p);
        }
        printf("\n");
    }
}

void dzero(double *v, int l);
void daxpy(double *restrict y, double a, const double *restrict x, int l);
double dnrm2(const double *x, int l);
void dscal(double *x, double a, int l);
double ddot(const double *x, const double *y, int l);
void dswap(double *x, double *y, int n);
void dcopy(double *x, double *y, int l);

void dzero(double *v, int l)
{
    memset(v, 0, sizeof(*v)*l);
}
void daxpy(double *restrict y, double a, const double *restrict x, int l)
{
        int m = l-3;
        int i;
        for (i = 0; i < m; i += 4)
        {
                y[i] += a * x[i];
                y[i+1] += a * x[i+1];
                y[i+2] += a * x[i+2];
                y[i+3] += a * x[i+3];
        }
        for ( ; i < l; ++i) /* clean-up loop */
                y[i] += a * x[i];
}
double dnrm2(const double *x, int l)
{
        double xx = ddot(x, x, l);
        return sqrt(xx);
}
void dscal(double *x, double a, int l)
{
        int m = l-4;
        int i;
        for (i = 0; i < m; i += 5){
                x[i] *= a;
                x[i+1] *= a;
                x[i+2] *= a;
                x[i+3] *= a;
                x[i+4] *= a;
        }

        for ( ; i < l; i++)        /* clean-up loop */
                x[i] *= a;
}
double ddot(const double *x, const double *y, int l)
{
        double s = 0;
        int m = l-4;
        int i;
        for (i = 0; i < m; i += 5)
                s += x[i] * y[i] + x[i+1] * y[i+1] + x[i+2] * y[i+2] +
                        x[i+3] * y[i+3] + x[i+4] * y[i+4];

        for ( ; i < l; i++)        /* clean-up loop */
                s += x[i] * y[i];

        return s;
}
void dswap(double *x, double *y, int n)
{
    for(int i=0; i<n; i++){
        double t = x[i];
        x[i] = y[i];
        y[i] = t;
    }
}
void dcopy(double *x, double *y, int l)
{
        memcpy(y, x, sizeof(*x)*(size_t)l);
}

int icmp(const void *a, const void *b)
{
    const Node *p = (const Node *)a, *q = (const Node *)b;
    if (p[0].index < q[0].index)
        return -1;
    else if (p[0].index > q[0].index)
        return 1;
    else if (p[1].index < q[1].index)
        return -1;
    else if (p[1].index > q[1].index)
        return 1;
    else
        return 0;
}

// Read the undiercted symmetric adjancency matrix
// into a csr format with duplicated edges (in place).
// Note that the order in each row is not guarenteed.
void read_undirected_adjmatrix(FILE *fin, CUTProblem *prob)
{
    int n, nnz;

    int buf_size = 120;
    char buf[buf_size];
    int is_mm_format = 0;
    int is_binary = 0;

    if(NULL == fgets(buf, buf_size, fin)) {
        fprintf(stderr, "Error: readline\n");
        exit(1);
    }
    if(buf[0] == '%') {
        is_mm_format = 1;
        if(NULL != strstr(buf, "pattern"))
            is_binary = 1;
    }
    fseek(fin, 0, SEEK_SET);
    
    int lineno = 1;
    // handle coments for matrix market
    while(NULL != fgets(buf, buf_size, fin)
            && buf[0] == '%')
        lineno++;

    if(is_mm_format) {
        int n_;
        int ret = sscanf(buf, " %d %d %d", &n, &n_, &nnz);
        if(ret != 3) {
            fprintf(stderr, "Error: matrix size\n");
            exit(1);
        }
        if(n != n_) {
            fprintf(stderr, "Error: not square matrix\n");
            exit(1);
        }
    } else {
        int ret = sscanf(buf, " %d %d", &n, &nnz);
        if(ret != 2) {
            fprintf(stderr, "Error: matrix size\n");
            exit(1);
        }
    }

    //printf("n=%d nnz=%d\n", n, nnz);

    prob->n = n;
    prob->C =    (Node**) Malloc(sizeof(Node*)*n);
    // Space usage: nnz*(double+int)*2 + nnz*int. 
    //              The nnz*int part is for construction, will free after use.
    // Time: O(nnz).
    // We store the COO format in ipart and jpart,
    // construct the CSR format for the upper triangle stored in COO,
    // allocate space for the transpose part,
    // and finally fill in the CSR format for the lower triangle.
    Node *pool = (Node*) Malloc(sizeof(Node)*(nnz*2+n));
    Node *jpart = pool+nnz+n;
    int *ipart = (int*) Malloc(sizeof(*ipart)*nnz);

    int *cap =      (int*) Calloc(n, sizeof(*cap));
    Node **rear   =    (Node**) Malloc (n*sizeof(Node*));
    Node end_node = {.index = -1, .value = 0};

    // read COO format
    int cur=0;
    while(1) {
        int i=0, j=0, t;
        double w=0;
        int ret, expect;

        if(NULL == fgets(buf, buf_size, fin)) {
            if(cur < nnz){
                fprintf(stderr, "ERROR less edges than specified\n");
                exit(1);
            }
            break;
        } else if (cur >= nnz) {
            fprintf(stderr, "ERROR more edges than specified\n");
            exit(1);
        }
        lineno++;
        
        if (is_binary) {
            ret = sscanf(buf, " %d %d", &i, &j);
            expect = 2;
            w = 1;
        } else {
            ret = sscanf(buf, " %d %d %lf", &i, &j, &w);
            expect = 3;
        }
        if (ret != expect) {
                fprintf(stderr, "ERROR reading file at %d\n", lineno);
                exit(1);
        }

        --i, --j;
        if (!(0<=i && i<n && 0<=j && j<n)) {
            fprintf(stderr, "ERROR at line %d\n", lineno);
            exit(1);
        } else if (i==j) {
            nnz--;
            continue;
        } else if (i>j) {
            t=i, i=j, j=t;
        }

        ipart[cur] = i;
        jpart[cur].index = j;
        jpart[cur].value = w;
        cur++, cap[i]++;
    }

    // construct CSR for the upper traingle
    rear[0] = pool;
    for (int i=1; i<n; i++) {
        rear[i] = rear[i-1] + cap[i-1];
    }
    for (cur=0; cur < nnz; cur++) {
        int i = ipart[cur], j = jpart[cur].index;
        *rear[i] = jpart[cur];
        rear[i]++;
        cap[j]++;
    }
    for (int i=0; i<n; i++) {
        cap[i]++; // add end_node
    }

    // expand CSR format
    Node* poolend  = pool+nnz*2+n;
    for (int i=n-1; i>0; i--) {
        Node* dst   = poolend - cap[i];
        Node* src   = rear[i-1];
        size_t len  = rear[i] - rear[i-1];

        memmove(dst, src, len*sizeof(Node));
        rear[i] += dst-src;
        poolend -= cap[i];
    }

    // fill in the lower triangle
    Node *front = pool;
    for (int i=0; i<n; i++) {
        for (Node *p=front; p<rear[i]; p++){
            int j = p->index;
            if (i>j) break;
            rear[j]->index = i;
            rear[j]->value = p->value;
            rear[j]++;
        }
        prob->C[i] = front;
        front += cap[i];
    }
    for (int i=0; i<n; i++) {
        *rear[i] = end_node;
    }
    prob->nnz = nnz;
    
    free(ipart);
    free(rear);
    free(cap);
}

void permutation(int *perm, int l)
{
        for(int i=0; i<l; i++){
                int j = (int)random() % (l-i);
                int t = perm[i];
                perm[i] = perm[j];
                perm[j] = t;
        }
}

void rand_unit(double *x, int k)
{
        for(int i=0; i<k; i++) {
                x[i] = drand48()*2-1;
        }
        double r = dnrm2(x, k);
        dscal(x, 1/r, k);
}

double square(double x)
{   return x*x;   }

double relu(double z)
{   
    double x = 1-z/4;
    return (x>=0) ? x : 0.; 
}

void eval_maxsat(SATProblem *prob, double **V, double **Z, int k, double *fval, int *rval, double *upper, Parameter *param)
{
    int m = prob->m, n = prob->n;
    int **index = prob->index;
    int *clen = prob->clause_len;
    int verbose = param->verbose;

    int *clause_sat = Calloc(m+1, sizeof(*clause_sat));
    double *r = Calloc(k, sizeof(*r));

    double ravg = 0;
    for(int trial=0; trial<param->n_trial; trial++){
        memset(clause_sat, 0, (m+1)*sizeof(*clause_sat));
        if(verbose) printf("v ");
        // rounding
        rand_unit(r, k);
        double wsign = ddot(r, V[0], k);
        //printf("v[0] = %f %f\n", V[0][0], V[0][1]);
        for(int j=1; j<=n; j++) {
            int bvar = (wsign * ddot(r, V[j], k) > 0) ? 1 : -1;
            int *p = index[j];
            for(; *p != 0; p++) {
                double sign = (*p > 0) ? 1. : -1;
                int i = abs(*p);
                if(sign * bvar > 0)
                    clause_sat[i] = sign*j;
            }
            if(verbose) printf("%d ", j * bvar);
        }
        if(verbose) printf("\n");

        if(verbose) printf("unsat clause ");
        *fval = *rval = *upper = 0;
        for(int i=1; i<=m; i++) {
            int div = (clen[i]+1)*clen[i]+2; // (len+1)^2-(len-1)
            double zi = square(dnrm2(Z[i], k));
            *fval += zi / div;
            if(isnan(*fval) || isinf(*fval)){
                fprintf(stderr, "zi %g i %d\n", zi, i);
                for(int kk=0; kk<k; kk++)
                    fprintf(stderr, "%f ", Z[i][kk]);
                fprintf(stderr, "\n");
                exit(1);
            }
            *rval += (clause_sat[i] != 0);
            *upper += 1-(zi-(clen[i]-1))/div;
            
            if(verbose && clause_sat[i] == 0) printf("%d ", i);
        }
        if(verbose) printf("\n");
        ravg += *rval;
    }
    ravg /= param->n_trial;
    *rval = ravg;

    free(clause_sat);
    free(r);
}

void do_maxsat(SATProblem *prob, Model *model, Parameter *param)
{
    int max_iter = param->max_iter;
    double eps = param->eps;

    int m = prob->m, n = prob->n;
    int k = (int)ceil(sqrt(2*n));
    if(param->k > 0){
        k = param->k;
    }else if(param->k < 0){
        k /= -param->k;
        if(k==0) k=1;
    }


    int **index = prob->index;
    int *clen = prob->clause_len;

    printf("MAXSAT m %d n %d k %d nnz %d eps %.4e max_iter %d\n", m, n, k, prob->nnz, eps, max_iter);

    double **V = Calloc(n+1, sizeof(*V));
    for(int j=0; j<=n; j++) {
        V[j] = Calloc(k, sizeof(**V));
    }
    double **Z = Calloc(m+1, sizeof(*Z));
    for(int i=0; i<=m; i++) {
        Z[i] = Calloc(k, sizeof(**Z));
    }
    double *g = Calloc(k, sizeof(*g));
    double *ci = Calloc(k, sizeof(*ci));
    double *r = Calloc(k, sizeof(*r));

    for(int j=0; j<=n; j++) {
        rand_unit(V[j], k);
    }

    // init Ci = \sum_j +-vj - (len(Ci)-1)w 
    for(int i=1; i<=m; i++) {
        dzero(g, k);
        daxpy(g, -1., V[0], k);
        daxpy(Z[i], 1., g, k);
    }
    for(int j=1; j<=n; j++) {
        int *p = index[j];
        for(; *p != 0; p++) {
            double sign = (*p > 0) ? 1. : -1;
            int i = abs(*p);
            daxpy(Z[i], sign, V[j], k);
        }
    }

    double fval=0, upper=0, fval_prev=0;
    int rval;
    eval_maxsat(prob, V, Z, k, &fval, &rval, &upper, param);
    printf("iter %3d  fval %.14e  delta %.4e  rval %d/%d  %.4e  upper %.4e  time %.14e\n",
            0, fval, 0., rval, m, rval*1./m, upper/m, 0.);

    fval_prev = fval;

    int64_t time_st = wall_clock_ns();
    int64_t time_eval = 0;
    int iter = 0;
    for(; iter < max_iter; iter++) {
        for(int j=1; j<=n; j++) {
            // Remove vj from Ci
            int *p;
            for(p = index[j]; *p != 0; p++) {
                double sign = (*p>0) ? 1. : -1.;
                int i = abs(*p);
                daxpy(Z[i], -sign, V[j], k);
            }
            // Do Mixing
            dzero(g, k);
            for(p = index[j]; *p; p++) {
                double sign = (*p>0) ? 1. : -1.;
                int i = abs(*p);
                daxpy(g, sign/((clen[i]+1)*clen[i]+2), Z[i], k);
            }
            double gnrm = dnrm2(g, k);
            if(gnrm >= MEPS){
                dscal(g, -1/gnrm, k);
                dcopy(g, V[j], k);
            }
            // Add vj back to Ci
            for(p = index[j]; *p != 0; p++) {
                double sign = (*p>0) ? 1. : -1.;
                int i = abs(*p);
                daxpy(Z[i], sign, V[j], k);
            }
        }
        int64_t time_now = wall_clock_ns();

        eval_maxsat(prob, V, Z, k, &fval, &rval, &upper, param);
        int64_t time_eval_ed = wall_clock_ns();
        time_eval += time_eval_ed-time_now;
        
        double delta = fval_prev - fval;
        printf("iter %3d  fval %.14e  delta %.4e  rval %d/%d  %.4e  upper %.4e  time %.14e\n",
                iter+1, fval, delta, rval, m, rval*1./m, upper/m, wall_time_diff(time_eval_ed-time_eval, time_st));
        fflush(stdout);
        if(delta < eps)
            break;
        fval_prev = fval;
    }

    model->n = n, model->k = k;
    model->V = V;
}

void do_maxcut(CUTProblem *prob, Model *model, Parameter *param)
{
    int iter;
    int max_iter = param->max_iter;
    int n = prob->n;
    int k = (int)ceil(sqrt(2*n));
    double eps = param->eps;
    if(param->k > 0){
        k = param->k;
    }else if(param->k < 0){
        k /= -param->k;
        if(k==0) k=1;
    }

    printf("MAXCUT n=%d k=%d nnz=%d eps %.4e max_iter %d\n", n, k, prob->nnz, eps, max_iter);

    double **V = (double**) Malloc (sizeof(*V)*n);
    V[0] = (double*) Calloc (n*k, sizeof(**V));

    for(int i=0; i<n; i++) {
            if(i!=0) V[i] = V[i-1] + k;
            rand_unit(V[i], k);
    }

    Node **C = prob->C;
    double *g = (double*) Malloc (sizeof(*g)*n);
    //int *perm = (int*) Malloc(sizeof(*perm)*n);
    //permutation(perm, n);

    // calculate function value
    double fval = 0;
    for(int i=0; i<n; i++) {
        for(const Node *p=C[i]; p->index != -1; p++) {
            fval += p->value * ddot(V[i], V[p->index], k);
            fval -= p->value;
        }
    }
    printf("iter %4d time %.16g\tfval %.22g\tdelta %.16g\n",
            0, 0.,
            fval/4, 0.);
    fflush(stdout);

    int64_t time_st = wall_clock_ns();

    for(iter=0; iter < max_iter; iter++) {
        double delta = 0;
        for(int i=0; i<n; i++) {
            //int i = ii; //perm[ii];
            dzero(g, k);
            for(const Node *p=C[i]; p->index != -1; p++) {
                daxpy(g, p->value, V[p->index], k);
            }
            double gnrm = dnrm2(g, k);
            if(gnrm > MEPS) {
                dscal(g, -1/gnrm, k);
            }
            delta += 2*gnrm*(1. - ddot(g, V[i], k));
            dcopy(g, V[i], k);
        }
        fval -= delta;
        int64_t time_now = wall_clock_ns();
        printf("iter %4d time %.16g\tfval %.22g\tdelta %.16g\n",
                iter+1,
                wall_time_diff(time_now, time_st),
                fval/4, delta/4);
        fflush(stdout);
        if(delta < eps)
            break;
    }

    model->n = n, model->k = k;
    model->V = V;
}

void get_clause_stat(SATProblem *prob, int *min_len, int *max_len, double *avg_len)
{
    int *len = prob->clause_len;
    *max_len = 0;
    *min_len = INT_MAX;
    *avg_len = 0.;
    for(int i=1; i<=prob->m; i++){
        if(*max_len < len[i])
            *max_len = len[i];
        if(*min_len > len[i])
            *min_len = len[i];
        *avg_len += len[i];
    }
    *avg_len /= prob->m;
}

void print_usage(char* prog_name, Parameter *param)
{
    printf( "%s [OPTIONS] INPUT: Mixing method for SDP\n", prog_name); 
    printf( "OPTIONS:\n");
    printf( "\t-s SOLVER: type of solver\n");
    printf( "\t           \"-s maxcut\" for maximum cut\n");
    printf( "\t           \"-s maxsat\" for maximum SAT (default)\n");
    printf( "\t-k RANK: rank of solution (default auto)\n");
    printf( "\t         use \"-k /2\" to divide the rank by 2\n");
    printf( "\t-e EPS: stopping threshold (default %.4e)\n", param->eps);
    printf( "\t-t MAX_ITER: maximum iteration (default %d)\n", param->max_iter);
    printf( "\t             use \"-t max\" for INT_MAX\n");
    printf( "\t-r N_TRIAL: number of trial in evaluation (default %d)\n", param->n_trial);
    printf( "\t-u: use unspeficied wcnf format\n");
    printf( "\t-v: verbose\n");
}

void get_parameter(int argc, char **argv, Parameter *param)
{
    Parameter _param = {
        .solver = MAXSAT,
		.k = -1,
        .eps = 1e-3,
        .max_iter = 1000,
        .is_unspeficied_wcnf = 0,
        .verbose = 0,
        .n_trial = 1,
        .fin = NULL,
        .fin_name = NULL,
    };

    if(argc <= 1){
        print_usage(argv[0], &_param);
        exit(0);
    }

    char **p = argv+1;
    int i;
    for(i=1; i<argc; i++, p++){
        if(!strcmp(*p, "-s")){
            if(i+1 >= argc) break;
            if(!strcmp(p[1], "maxcut")){
                _param.solver = MAXCUT;
            }else if(!strcmp(p[1], "maxsat")){
                _param.solver = MAXSAT;
            }else {
                int ret = sscanf(p[1], "%d", &_param.solver);
                if(ret != 1 || !(_param.solver >=0 && _param.solver <= 1)) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-k")){
            if(i+1  >= argc) break;
            int ret = sscanf(p[1], "/%d", &_param.k);
            if(ret==1){
                _param.k *= -1;
            }else{
                ret = sscanf(p[1], "%d", &_param.k);
                if(ret != 1 || _param.k <= 0) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-e")){
            if(i+1 >= argc) break; 
            int ret = sscanf(p[1], "%lf", &_param.eps);
            if(ret != 1) break; 
            i++, p++;
        }else if(!strcmp(*p, "-t")){
            if(i+1 >= argc) break; 
            if(!strcmp(p[1], "max")){
                _param.max_iter = INT_MAX;
            }else{
                int ret = sscanf(p[1], "%d", &_param.max_iter);
                if(ret != 1) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-r")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%d", &_param.n_trial);
            if(ret != 1) break; 
            if(_param.n_trial < 1)
                _param.n_trial = 1;
            i++, p++;
        }else if(!strcmp(*p, "-u")){
            _param.is_unspeficied_wcnf = 1;
        }else if(!strcmp(*p, "-v")){
            _param.verbose = 1;
        }else if(i+1 == argc){
            _param.fin = fopen(*p, "r");
            if(!_param.fin){
                fprintf(stderr, "%s\n", strerror(errno));
                exit(1);
            }
            _param.fin_name = strdup(*p);
        }else{
            printf("Error: no such parameter\n");
            break;
        }
    }
    if(i != argc || !_param.fin){
        print_usage(argv[0], &_param);
        exit(0);
    }
    *param = _param;
}

int main(int argc, char **argv)
{
    Parameter param;
    get_parameter(argc, argv, &param);
    srand48(0);

    Model model;
    if(param.solver == MAXCUT){
        printf("Solving maximum cut. Reading %s\n", param.fin_name);
        CUTProblem cut_prob;
        read_undirected_adjmatrix(param.fin, &cut_prob);
        do_maxcut(&cut_prob, &model, &param);
    }else if(param.solver == MAXSAT){
        SATProblem sat_prob;
        printf("Solving maximum SAT. Reading %s\n", param.fin_name);
        read_sat_problem(param.fin, &sat_prob, param.is_unspeficied_wcnf);
        if(param.verbose)
            print_struct(&sat_prob);

        int min_len, max_len;
        double avg_len;
        get_clause_stat(&sat_prob, &min_len, &max_len, &avg_len);
        printf("clauses len : min %d max %d avg %.4f\n", min_len, max_len, avg_len);

        do_maxsat(&sat_prob, &model, &param);
    }

    FILE *fout = fopen(strcat(param.fin_name, ".sol"), "w");
    for (int i=0; i<model.n; i++) {
        for (int j=0; j<model.k; j++)
            fprintf(fout, "%.22g ", model.V[i][j]);
        fprintf(fout, "\n");
    }

    return 0;
}
