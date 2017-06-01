# Usage

```
./mixing [OPTIONS] DATA_FILE
OPTIONS:
	-s SOLVER: type of solver
	           "-s maxcut" for maximum cut
	           "-s maxsat" for maximum SAT (default)
	-k RANK: rank of solution (default auto)
	         use "-k /2" to divide the rank by 2
	-e EPS: stopping threshold (default 1.0000e-03)
	-t MAX_ITER: maximum iteration (default 1000)
	             use "-t max" for INT_MAX
	-r N_TRIAL: number of trial in evaluation (default 1)
	-u: use unspeficied wcnf format
	-v: verbose

To complile the file, please use
	$ make
```

For any questions and comments, please send your email to
poweiw@cs.cmu.edu

