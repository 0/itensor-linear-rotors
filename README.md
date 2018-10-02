# itensor-linear-rotors

An [ITensor](http://itensor.org/) [`SiteSet`](http://itensor.org/docs.cgi?page=classes/siteset) for linear rotors.


## Examples

* `bin/dipoles_dmrg 1.0 4 3 data/sample_sweep_table 2e-4 5 20`

---

* `bin/dipoles_dmrg_write 1.0 4 3 data/sample_sweep_table 2e-4 5 20 sites H mps`
* `bin/dipoles_dmrg_read 4 3 sites H mps`

---

* `bin/dipoles_dmrg_probabilities 1.0 3 3 data/sample_sweep_table 2e-4 5 20`
* `bin/dipoles_dmrg_sampling 1.0 3 3 data/sample_sweep_table 2e-4 5 20 100`

---

* `bin/dipoles_dmrg_field 1.0 4 0.02 3 data/sample_sweep_table 2e-4 5 20`
* `bin/dipoles_dmrg_nonlinear data/geom_tetrahedron 1.5 4 2 data/sample_sweep_table 2e-4 5 20`


## Publications

The following publications contain data created using this package:

* Dmitri Iouchtchenko and Pierre-Nicholas Roy. **Ground states of linear rotor chains via the density matrix renormalization group.** The Journal of Chemical Physics 148, 134115 (2018). [doi:10.1063/1.5024403](https://aip.scitation.org/doi/abs/10.1063/1.5024403), [arXiv:1805.05420](https://arxiv.org/abs/1805.05420).


## References

The sampling algorithm is from [Ferris and Vidal (2012) "Perfect sampling with unitary tensor networks"](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.165146) ([preprint](https://arxiv.org/abs/1201.3974)).


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
