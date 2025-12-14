# ZDF - Run Instructions

Before a fresh filtering run, you'll need to be clear on:
* `N` The length of your filter. `N ≤ 65535`.
* `{d_1, d_2, ..., d_m}` The list of derivatives to filter for. `d_i < 16`.
* `q` The degree of non-minimality of the filter. `q < 16`.
* `K` The shape parameter $\kappa$, usually zero. `K < 16`.
* `M` The lowest $\mu$ shape parameter to filter with. `M < 16`. 
* `U` The highest-plus-one $\mu$ shape parameter to filter with. `U < 16`. 

The filtering output, when run against a data file of length `M`, is an `M x ((U-M).m)` array organized as follows: 
$$
\left[\begin{array}{ccccccc} d_{1}^{\mu_1} & d_{1}^{\mu_2} & ... & d_{1}^{\mu_{U-M}} & d_{2}^{\mu_1} & ... & d_{m}^{\mu_{U-M}} \end{array}\right]
$$

1. Retrieve the encoding of the derivatives by concatenating `d_1, d_2, ..., d_m`, in that order, as an argument to `./zdf -d `. For example: ```./zdf -d 0,1,2```

1. Armed with the resulting integer `D`, say `7`, similarly retrieve the encoding of the whole ZDF by using flag `./zdf -e N,D,q,K,M,U`. For example: ```./zdf -e 1024,7,2,0,1,2```

1. Use the resulting integer `T` to label your initialization file. That should be a binary file with exactly `N` floats written to it, labeled `${T}.zdft` and placed in your data directory. See [this example](../tst/TestData.py). Double check that the `REPO` and `T` constants in [zdf.cpp](../src/zdf.cpp) match your data directory and encoding. 

1. If using the executable to filter observations from a file: That file should also be generated as `${T}.zdfi`. Ensure that the update line in [zdf.cpp](../src/zdf.cpp) reads `zdf.update<L>(REPO);` where `L` is at least the number of entries in your `.zdfi` file. The filtered array is output to `${T}.zdfo`.

1. Build the executable, similarly to the BPPR [build instructions](https://github.com/gcbeck/bppr/blob/master/doc/build.md)

1. Run the executable against your file(s). For example, for a single MKL thread and to specify that the initialization file should be overwritten with any new update observations: 
```
./zdf -t 1 -w
```

