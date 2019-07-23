# itensor-linear-rotors

An [ITensor](http://itensor.org/) [`SiteSet`](http://itensor.org/docs.cgi?page=classes/siteset) for linear rotors.


## Examples

* `bin/dipoles_dmrg -R 1.0 -N 4 --l-max 3 --sweep-table data/sample_sweep_table --dH2-goal 2e-4 --sweeps-min 5 --sweeps-max 20`

---

* `bin/dipoles_dmrg_write -R 1.0 -N 4 --l-max 3 --sweep-table data/sample_sweep_table --dH2-goal 2e-4 --sweeps-min 5 --sweeps-max 20 --sites-path data/sites --H-path data/H --mps-path data/mps`
* `bin/dipoles_dmrg_read --sites-path data/sites --H-path data/H --mps-path data/mps`

---

* `bin/dipoles_dmrg_probabilities -R 1.0 -N 3 --l-max 3 --sweep-table data/sample_sweep_table --dH2-goal 2e-4 --sweeps-min 5 --sweeps-max 20`
* `bin/dipoles_dmrg_sampling -R 1.0 -N 3 --l-max 3 --sweep-table data/sample_sweep_table --dH2-goal 2e-4 --sweeps-min 5 --sweeps-max 20 --num-samples 100`

---

* `bin/dipoles_dmrg_field -R 1.0 --field 0.02 -N 4 --l-max 3 --sweep-table data/sample_sweep_table --dH2-goal 2e-4 --sweeps-min 5 --sweeps-max 20`
* `bin/dipoles_dmrg_nonlinear --geom data/geom_tetrahedron -R 1.5 -N 4 --l-max 2 --sweep-table data/sample_sweep_table --dH2-goal 2e-4 --sweeps-min 5 --sweeps-max 20`


## Publications

The following publications contain data created using this package:

* Dmitri Iouchtchenko and Pierre-Nicholas Roy. **Ground states of linear rotor chains via the density matrix renormalization group.** The Journal of Chemical Physics 148, 134115 (2018). [doi:10.1063/1.5024403](https://aip.scitation.org/doi/abs/10.1063/1.5024403), [arXiv:1805.05420](https://arxiv.org/abs/1805.05420).


## References

The sampling algorithm is from [Ferris and Vidal (2012) "Perfect sampling with unitary tensor networks"](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.165146) ([preprint](https://arxiv.org/abs/1201.3974)).


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
