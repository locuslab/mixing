# Mixing
The Mixing method for maximum cut (MAXCUT) and maximum satisfiability (MAXSAT) problem.

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
```

To compile the file, please use
```
	$ make
```

# More Info
This repository is by [Po-Wei Wang](http://powei.tw),
[Wei-Cheng Chang](https://octoberchang.github.io/),
and [J. Zico Kolter](http://zicokolter.com)
and contains the source code to
reproduce the experiments in our paper
[The Mixing method: coordinate descent for low-rank semidefinite programming](http://arxiv.org/abs/1706.00476).
If you find this repository helpful in your publications, please consider citing our paper.
```
@article{wang2017mixing,
	title = {The Mixing method: coordinate descent for low-rank semidefinite programming},
	author = {Po-Wei Wang and Wei-Cheng Chang and J. Zico Kolter},
	journal = {arXiv preprint arXiv:1706.00476},
	year = {2017}
}
```

For any questions and comments, please send your email to
[poweiw@cs.cmu.edu](mailto:poweiw@cs.cmu.edu)

