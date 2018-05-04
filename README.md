# Counting Linear Extensions in Practice

The code used for the experiments in the IJCAI 2018 paper A Scalable Scheme for Counting Linear Extensions and AAAI 2018 paper Counting Linear Extensions in Practice: MCMC versus Exponential Monte Carlo

## COMPILING

The Makefile in the root can be used to compile all solutions and generate all
instances:

```$ make```

The running make in the subdirectories can be used to compile parts
selectively.

## DEPENDENCIES

Instance generation requires Python 3. The solutions require a version of the
g++ compiler that supports C++11. In addition, the LEcount solutions (Exact
Dynamic Programming and Adaptive Relaxation Monte Carlo) require Boost and GMP.

## USAGE

The instances are generated to instances/*.txt as adjacency matrices. The
draw.py script can be used to visualize them to instances/img/ (requires
GraphViz).

Each solution program in solutions/ takes the name of the adjacency matrix file
as the first command line argument, possibly followed by solution-specific
options, and outputs the estimate for the number of linear extensions. Some
solutions also write extra information to the standard error stream. All
solutions support partial orders with at most 512 elements.

### IJCAI 2018: A Scalable Scheme for Counting Linear Extensions

Relaxation Tootsie Pop:

```$ solutions/relaxation-tpa/count INSTANCE```

Trivial Relaxation Tootsie Pop:

```$ solutions/relaxation-tpa/trivial INSTANCE```

Telescopic Product:

```$ solutions/telescopic-product/decomposition INSTANCE GibbsLinextSampler```

Extension Tootsie Pop:

```$ solutions/extension-tpa/count INSTANCE```

Adaptive Relaxation Monte Carlo:

```$ solutions/lecount/lecount INSTANCE --algorithm=armc```

SAT Encodings #1 and #2:

These output a DIMACS CNF encoding for the instance, which can then be used as input for a SAT model counter (D4, sharpSAT and ApproxMC2 in the paper) to count the linear extensions.

```$ solutions/sat/encoding.py INSTANCE 1```
```$ solutions/sat/encoding.py INSTANCE 2```

### AAAI 2018: Counting Linear Extensions in Practice: MCMC versus Exponential Monte Carlo

Telescopic Product:

```$ solutions/telescopic-product/basic INSTANCE SwapLinextSampler```

Decomposition Telescopic Product:

```$ solutions/telescopic-product/decomposition INSTANCE SwapLinextSampler```

Decomposition Telescopic Product using the Gibbs sampler:

```$ solutions/telescopic-product/decomposition INSTANCE GibbsLinextSampler```

Tootsie Pop:

```$ solutions/tpa/count INSTANCE```

Exact Dynamic Programming:

```$ solutions/lecount/lecount INSTANCE --algorithm=dp```

Adaptive Relaxation Monte Carlo:

```$ solutions/lecount/lecount INSTANCE --algorithm=armc```

Variable Elimination via Inclusion-Exclusion (exact):

```$ solutions/lecount/lecount INSTANCE --algorithm=veie```
