# itensor-linear-rotors

An [ITensor](http://itensor.org/) [`SiteSet`](http://itensor.org/docs.cgi?page=classes/siteset) for linear rotors.


## Example workflow

### Generate sites and Hamiltonian MPO

#### Full symmetry

1. `bin/gen_sites -N 4 --l-max 2 --sites-out-path data/sites`
1.
   * `bin/gen_H -R 1.0 --mpo-cutoff 1e-20 --sites-in-path data/sites --H-out-path data/H`
   * `bin/gen_H -g 1.0 --sociability 1 --mpo-cutoff 1e-20 --sites-in-path data/sites --H-out-path data/H`
   * `bin/gen_H --pbc -R 1.0 --mpo-cutoff 1e-20 --sites-in-path data/sites --H-out-path data/H`

#### Reduced symmetry

1. `bin/gen_sites -N 4 --l-max 2 --m-sym 0 --sites-out-path data/sites`
1. `bin/gen_H --geom-in-path data/geom_tetrahedron -R 1.0 --mpo-cutoff 1e-20 --sites-in-path data/sites --H-out-path data/H`

#### No symmetry

1. `bin/gen_sites -N 4 --l-max 2 --lp-sym 0 --sites-out-path data/sites`
1. `bin/gen_H --field 0.02 --field-linear -R 1.0 --mpo-cutoff 1e-20 --sites-in-path data/sites --H-out-path data/H`

### Sweep MPS with DMRG

1. `bin/dmrg_sweep --sweep-table data/sample_sweep_table --num-sweeps 5 --sites-in-path data/sites --H-in-path data/H --mps-out-path data/mps05`
1. `bin/dmrg_sweep --sweep-table data/sample_sweep_table --num-sweeps 5 --skip-sweeps 5 --sites-in-path data/sites --H-in-path data/H --mps-in-path data/mps05 --mps-out-path data/mps10`

### Perform analysis

* `bin/analyze --sites-in-path data/sites --mps-in-path data/mps10`
* `bin/sample --num-samples 100 --sites-in-path data/sites --mps-in-path data/mps10`
* `bin/dump_state --sites-in-path data/sites --mps-in-path data/mps10`


## Publications

The following publications contain data created using this package:

* Dmitri Iouchtchenko and Pierre-Nicholas Roy. **Ground states of linear rotor chains via the density matrix renormalization group.** The Journal of Chemical Physics 148, 134115 (2018). [doi:10.1063/1.5024403](https://aip.scitation.org/doi/abs/10.1063/1.5024403), [arXiv:1805.05420](https://arxiv.org/abs/1805.05420).


## References

The sampling algorithm is from [Ferris and Vidal (2012) "Perfect sampling with unitary tensor networks"](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.165146) ([preprint](https://arxiv.org/abs/1201.3974)).


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
