# Base-2 Differential Privacy Reference Implementations
Reference implementations for base-2 differential privacy.

Author: Christina Ilvento


### Papers and other reference materials
* *Base-2 DP and Base-2 Exponential Mechanism*, [Ilvento '19](https://arxiv.org/abs/1912.04222)
* [Abstract](https://drive.google.com/file/d/1OytgB24d1n-xPIWrrKCsVQQdS7rV3tjn/view?usp=sharing) - *Implementing Sparse Vector with Base-2 DP*
* [Abstract](https://drive.google.com/file/d/1okHAkjNENiS2WfSKdkUo8B29yE8-Qfof/view?usp=sharing) - *Implementing differentially private integer partitions with Base-2 DP* an extension of [Blocki, Datta and Bonneau '16](http://www.jbonneau.com/doc/BDB16-NDSS-pw_list_differential_privacy.pdf).

### [Documentation](https://cilvento.github.io/b2dp/doc/b2dp/index.html)
* [Example notebook](https://bit.ly/3nviEiD) for sparse-vector technique issues.

### Benchmarks
[`Criterion`](https://github.com/bheisler/criterion.rs) benchmarks can be found in `b2dp/src/benchmarks/`. To run the benchmarks: 
```bash
cd rust/b2dp
cargo bench
```
The benchmark results can then be found in the `rust/b2dp/target/Criterion` directory.

### Python Version
An older Python version of the base-2 exponential mechanism along with attack demonstrations is archived in the [b2_exponential_mechanism](https://github.com/cilvento/b2_exponential_mechanism) repository. Future development is planned in Rust.

### Dependencies
* [GMP](https://gmplib.org/manual/Installing-GMP.html), [MPFR](https://www.mpfr.org/mpfr-current/mpfr.html)
* [rug](https://docs.rs/rug/1.8.0/rug/), [gmp-mpfr-sys](https://docs.rs/gmp-mpfr-sys/1.2.2/gmp_mpfr_sys/index.html)
* (for benchmarks) [Criterion](https://github.com/bheisler/criterion.rs)
