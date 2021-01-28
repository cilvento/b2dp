initSidebarItems({"fn":[["adjust_eta","Adjust the given eta to account for granularity 1/gamma_inv. "],["conditional_lazy_threshold","Determines whether the discretized Laplace exceeds the given  threshold conditioned on it exceeding the given conditional threshold. ## Arguments   * `eta`: the privacy parameter   * `arithmeticconfig`: ArithmeticConfig with sufficient precision   * `gamma`: granularity parameter   * `threshold`: the threshold value    * `cond_threshold`: the conditional threshold value already exceeded      (Must be smaller than `threshold`)   * `rng`: randomness source   * `optimize`: whether to optimize sampling, exacerbates timing channels"],["get_sum","Returns the sum:  (1-base)*\\sum_{k=start}^{end}base^{|k|}"],["is_multiple_of","Returns whether a is an integer multiple of b."],["lazy_threshold","Determines whether the discretized Laplace exceeds the given threshold.  ## Arguments   * `eta`: the privacy parameter   * `arithmeticconfig`: ArithmeticConfig with sufficient precision   * `gamma`: granularity parameter   * `threshold`: the threshold value    * `rng`: randomness source   * `optimize`: whether to optimize sampling, exacerbates timing channels"],["sample_within_bounds","Samples from the Discrete Laplace mechanism at granularity `gamma` within the provided bounds `wmin` and `wmax` where the values `wmin` and `wmax` are sampled with probability equal to the sum of the probabilities of all elements greater than or less than the bounds, respectively.  ## Arguments   * `eta`: privacy parameter  * `gamma`: granularity parameter, must be reciprocal of an integer.   * `wmin`: the minimum bound, must be an integer multiple of `gamma`.  * `wmax`: the maximum bound, must be an integer multiple of `gamma`.  * `arithmeticconfig`: ArithmeticConfig with sufficient precision for      sampling.   * `rng`: randomness source  * `optimize`: whether to optimize sampling. "]]});