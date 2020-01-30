# itensor-linear-rotors

An [ITensor](http://itensor.org/) [`SiteSet`](http://itensor.org/docs.cgi?page=classes/siteset) for linear rotors.


## Configuration

The `Makefile` assumes that you have the following standard environment variables set:

* `CPATH` should contain the path to the ITensor root directory.
* `LIBRARY_PATH` should contain the path to ITensor's `lib` directory.

Additionally, the options chosen in ITensor's `options.mk` file will determine the values needed in the following environment variables:

* `ITENSOR_CXX` should be set to the compiler used to build ITensor (likely `g++` or `clang++`).
* `ITENSOR_CXXFLAGS` should be set to the list of compiler flags needed by ITensor (for example, `-std=c++17 -DPLATFORM_openblas`)
* `ITENSOR_LDLIBS` should be set to the list of libraries to link against (possibly `-litensor -lpthread -lopenblas`).


### Linux, GCC, OpenBLAS

Sample configuration for Linux with [GCC](https://gcc.gnu.org/) and [OpenBLAS](https://www.openblas.net/):

```
export ITENSOR_PATH=...
export CPATH="${ITENSOR_PATH}:${CPATH}"
export LIBRARY_PATH="${ITENSOR_PATH}/lib:${LIBRARY_PATH}"
export ITENSOR_CXX=g++
export ITENSOR_CXXFLAGS='-std=c++17 -DPLATFORM_openblas -fpermissive -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE'
export ITENSOR_LDLIBS='-litensor -lpthread -lopenblas'
```

### macOS, Clang

Sample configuration for macOS with [Clang](https://clang.llvm.org/):

```
export ITENSOR_PATH=...
export CPATH="${ITENSOR_PATH}:${CPATH}"
export LIBRARY_PATH="${ITENSOR_PATH}/lib:${LIBRARY_PATH}"
export ITENSOR_CXX=clang++
export ITENSOR_CXXFLAGS='-std=c++17 -DPLATFORM_macos'
export ITENSOR_LDLIBS='-litensor -framework Accelerate'
```


## Example workflow

### Generate sites and Hamiltonian MPO

#### Full symmetry

1. `bin/gen_sites -N 4 --l-max 2 --sites-out workspace/sites`
1.
   * `bin/gen_ham -R 1.0 --mpo-cutoff 1e-20 --sites workspace/sites --ham-out workspace/ham`
   * `bin/gen_ham --anisotropy 0.0 -g -1.0 --sociability 1 --mpo-cutoff 1e-20 --sites workspace/sites --ham-out workspace/ham`
   * `bin/gen_ham --pbc -R 1.0 --mpo-cutoff 1e-20 --sites workspace/sites --ham-out workspace/ham`

#### Reduced symmetry

1. `bin/gen_sites -N 4 --l-max 2 --m-sym 0 --sites-out workspace/sites`
1. `bin/gen_ham --geom data/geom_tetrahedron -R 1.0 --mpo-cutoff 1e-20 --sites workspace/sites --ham-out workspace/ham`

#### No symmetry

1. `bin/gen_sites -N 4 --l-max 2 --lp-sym 0 --sites-out workspace/sites`
1. `bin/gen_ham --field 0.02 --field-linear -R 1.0 --mpo-cutoff 1e-20 --sites workspace/sites --ham-out workspace/ham`

### Sweep MPS with DMRG

1. `bin/dmrg_sweep --sweep-table data/sample_sweep_table_fast --num-sweeps 5 --sites workspace/sites --ham workspace/ham --mps-out workspace/mps05`
1. `bin/dmrg_sweep --sweep-table data/sample_sweep_table_fast --mps-cutoff 1e-14 --num-sweeps 5 --first-sweep 6 --sites workspace/sites --ham workspace/ham --mps-in workspace/mps05 --mps-out workspace/mps10`
1. `bin/dmrg_sweep --sweep-table data/sample_sweep_table_slow --num-sweeps 10 --sites workspace/sites --ham workspace/ham --ortho-mps workspace/mps10 --mps-out workspace/mps_exc`

### Perform analysis

* `bin/analyze --sites workspace/sites --mps workspace/mps10`
* `bin/sample --num-samples 100 --sites workspace/sites --mps workspace/mps10`
* `bin/dump_state --sites workspace/sites --mps workspace/mps10`

### Miscellaneous operations

* `bin/embiggen --sites1 workspace/sites --mps1-in workspace/mps10 --sites2 workspace/sites_big --mps2-out workspace/mps10_big`


## Publications

The following publications contain data created using this package:

* Dmitri Iouchtchenko and Pierre-Nicholas Roy. **Ground states of linear rotor chains via the density matrix renormalization group.** The Journal of Chemical Physics 148, 134115 (2018). [doi:10.1063/1.5024403](https://aip.scitation.org/doi/abs/10.1063/1.5024403), [arXiv:1805.05420](https://arxiv.org/abs/1805.05420).


## References

The sampling algorithm is from [Ferris and Vidal (2012) "Perfect sampling with unitary tensor networks"](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.165146) ([preprint](https://arxiv.org/abs/1201.3974)).


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
