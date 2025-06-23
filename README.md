# MultiSpinCoding

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ArrogantGao.github.io/MultiSpinCoding.jl/stable/) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ArrogantGao.github.io/MultiSpinCoding.jl/dev/) -->
[![Build Status](https://github.com/ArrogantGao/MultiSpinCoding.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/MultiSpinCoding.jl/actions/workflows/CI.yml?query=branch%3Amain)
<!-- [![Coverage](https://codecov.io/gh/ArrogantGao/MultiSpinCoding.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ArrogantGao/MultiSpinCoding.jl) -->


An implementation of multi-spin coding for simulating the spin-glass model, according to "Optimised simulated annealing for Ising spin glasses".
The one implemented here is the method named `an_ms_r1_nf` in the paper, which requires $s = \pm 1$, $J_{ij} = \pm 1$ and $h = 0$.

## Usage

First install the package by typing the following in the Julia REPL:
```julia
] add https://github.com/ArrogantGao/MultiSpinCoding.jl
```

A simple example:
```julia
using MultiSpinCoding

# load the graph and the coupling matrix
g, J = load_instance("example/square16x16.txt")

# generate the model, n = 256 is the number of replicas
model = SpinGlass(g, J, 256)

# generate the scheduler, 2000 sweeps, beta0 = 0.1, beta1 = 3.0
scheduler = Scheduler(2000, 0.1, 3.0)

# run the simulated annealing
sa!(model, scheduler)

# calculate the energies
energies = cal_energies(model)
```

One can also use the script `example/square16x16.jl` to run the example, which will generate a histogram of the final energy of the replicas. The result is shown below:

![](example/square16x16.png)

## Benchmark

A simple benchmark is shown below:

```julia
julia> using MultiSpinCoding, BenchmarkTools

julia> g, J = load_instance("example/square16x16.txt");

julia> scheduler = Scheduler(2000, 0.1, 3.0);

# 64 replicas stored in a single UInt64
julia> model = SpinGlass(g, J, 64);

julia> @benchmark sa!($model, $scheduler)
BenchmarkTools.Trial: 764 samples with 1 evaluation per sample.
 Range (min … max):  6.038 ms …   6.982 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     6.236 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   6.506 ms ± 345.907 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

           ▂█▂                                            ▇█   
  ▃▃▃▂▂▂▂▄████▇▆▃▃▂▂▂▁▂▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▁▁▂▁▁▁▁▃▃▆▇██▇ ▃
  6.04 ms         Histogram: frequency by time        6.92 ms <

 Memory estimate: 0 bytes, allocs estimate: 0.

# 1024 replicas stored in a LongLongUInt{16}
julia> model = SpinGlass(g, J, 1024);

julia> @benchmark sa!($model, $scheduler)
BenchmarkTools.Trial: 368 samples with 1 evaluation per sample.
 Range (min … max):  12.192 ms …  14.028 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     13.847 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   13.585 ms ± 522.950 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

                                                       ▇▁█▅▆▆▁  
  ▂▂▁▁▁▁▁▂▁▁▁▁▃█▆▇▆▄▄▄▁▂▁▁▁▂▁▁▁▁▃▁▃▃▁▁▂▃▁▂▁▃▁▁▁▁▂▃▄▃▂▃████████ ▃
  12.2 ms         Histogram: frequency by time           14 ms <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> 6.038 / 2000 / 64 * 1e6
47.171875 # 47.171875 ns per round of update    

julia> 12.192 / 2000 / 1024 * 1e6
5.953125 # 5.953125 ns per round of update
```