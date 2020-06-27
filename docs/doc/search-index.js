var searchIndex={};
searchIndex["b2dp"] = {"doc":"Base-2 Differential Privacy Crate Implements the…","i":[[0,"utilities","b2dp","Base-2 Differential Privacy Utilities",null,null],[0,"params","b2dp::utilities","Implements the base-2 privacy parameter `eta`, which takes…",null,null],[3,"Eta","b2dp::utilities::params","Privacy parameter of the form `Eta = -z * log_2(x/2^y)`…",null,null],[12,"x","","",0,null],[12,"y","","",0,null],[12,"z","","",0,null],[11,"new","","Creates `Eta` privacy parameter from the given `x`, `y`…",0,[[["u32"]],[["eta"],["result",["eta","str"]],["str"]]]],[11,"from_epsilon","","Creates `Eta` privacy parameter approximating the base `e`…",0,[[["f64"]],[["eta"],["result",["eta","str"]],["str"]]]],[11,"check","","",0,[[["self"]],[["result",["str"]],["str"]]]],[11,"get_approximate_epsilon","","Returns the approximate corresponding epsilon Note that…",0,[[["self"]],["f64"]]],[11,"get_base","","Get the base `2^(-eta)` ## Arguments * precision: the…",0,[[["u32"],["self"]],[["float"],["str"],["result",["float","str"]]]]],[0,"exactarithmetic","b2dp::utilities","Implements methods requiring exact arithmetic, and…",null,null],[3,"ArithmeticConfig","b2dp::utilities::exactarithmetic","The exact arithmetic configuration. Includes the precision…",null,null],[12,"precision","","The required precision (computed based on other parameters)",1,null],[12,"inexact_arithmetic","","Whether an inexact operation has been performed in the…",1,null],[5,"randomized_round","","Randomized Rounding ## Arguments * `x`: the value to round…",null,[[["arithmeticconfig"],["f64"],["threadrandgen"]],["i64"]]],[5,"normalized_sample","","Normalized Weighted Sampling Returns the index of the…",null,[[["threadrandgen"],["arithmeticconfig"],["vec"],["bool"]],[["result",["usize","str"]],["str"],["usize"]]]],[11,"basic","","A basic arithmetic_config with default precision",1,[[],[["str"],["result",["arithmeticconfig","str"]],["arithmeticconfig"]]]],[11,"for_exponential","","Initialize an ArithmeticConfig for the base-2 exponential…",1,[[["i64"],["u32"],["eta"],["bool"]],[["str"],["result",["arithmeticconfig","str"]],["arithmeticconfig"]]]],[11,"increase_precision","","Increase the precision by `increment`. Returns an error…",1,[[["self"],["u32"]],[["result",["str"]],["str"]]]],[11,"check_mpfr_flags","","Check the current state of the flags",1,[[],[["result",["str"]],["str"]]]],[11,"set_retries","","Set the minimum number of retries for timing channel…",1,[[["self"],["u32"]],[["result",["str"]],["str"]]]],[11,"is_valid","","Returns whether the config is currently in a valid state",1,[[["self"]],["bool"]]],[11,"invalidate","","Invalidates the config",1,[[["self"]]]],[11,"enter_exact_scope","","Enter exact arithmetic scope. This method clears `mpfr`…",1,[[["self"]],[["result",["str"]],["str"]]]],[11,"exit_exact_scope","","Exit the exact arithmetic scope. Must be called after any…",1,[[["self"]],[["result",["str"]],["str"]]]],[0,"randomness","b2dp::utilities","Implements `ThreadRandGen` wrapper(s) for randomness…",null,null],[3,"GeneratorOpenSSL","b2dp::utilities::randomness","Random Number Generator using OpenSSL for randomness",null,null],[0,"discretesampling","b2dp::utilities","Implements discrete sampling methods in base-2 including…",null,null],[5,"is_multiple_of","b2dp::utilities::discretesampling","Returns whether a is an integer multiple of b",null,[[["float"]],["bool"]]],[5,"sample_within_bounds","","Samples from the Discrete Laplace mechanism at granularity…",null,[[["eta"],["threadrandgen"],["float"],["arithmeticconfig"],["copy"],["bool"]],[["float"],["str"],["result",["float","str"]]]]],[5,"adjust_eta","","Adjust the given eta to account for granularity 1/gamma_inv",null,[[["arithmeticconfig"],["eta"],["float"]],[["eta"],["result",["eta","str"]],["str"]]]],[5,"noisy_threshold","","Determines whether the discretized Laplace exceeds the…",null,[[["eta"],["threadrandgen"],["arithmeticconfig"],["float"],["bool"]],[["float"],["str"],["result",["float","str"]]]]],[5,"get_sum","","Returns the sum:…",null,[[["arithmeticconfig"],["float"],["float"]],[["float"],["str"],["result",["float","str"]]]]],[0,"bounds","b2dp::utilities","Compute partition bounds for the integer partition…",null,null],[3,"PartitionBoundOptions","b2dp::utilities::bounds","",null,null],[12,"exp_options","","Options to pass for any executions of the exponential…",2,null],[12,"sparsity_control","","Optional budget to use for sparsity control.",2,null],[12,"sparsity_padding","","Padding coefficient to use for sparsity control…",2,null],[12,"estimate_padding","","Padding for noisy estimates (multiplied by 1/epsilon)",2,null],[3,"PartitionBound","","Bounds for an integer partition exponential mechanism…",null,null],[12,"upper","","The upper bounds (terminates with 0)",3,null],[12,"lower","","The lower bounds (terminates with 0)",3,null],[12,"count","","Total count (bound)",3,null],[12,"cells","","Number of cells (bound)",3,null],[12,"sparsity_control","","Whether sparsity control is enforced",3,null],[12,"noisy_estimates","","The noisy estimates used, if any",3,null],[11,"check","","Checks that the partition bounds are valid: * Upper bound…",3,[[["self"]],[["result",["str"]],["str"]]]],[11,"new","","Get partition bounds for the given total count based on…",3,[[["usize"]],[["partitionbound"],["str"],["result",["partitionbound","str"]]]]],[11,"with_cells","","Get partition bounds for the given by the naive formula…",3,[[["usize"]],[["partitionbound"],["str"],["result",["partitionbound","str"]]]]],[11,"from_dist","","Returns partition bounds of distance at most `d` from the…",3,[[["vec"],["usize"]],[["partitionbound"],["str"],["result",["partitionbound","str"]]]]],[11,"with_reference","","Get bounds based on a reference partition using budget to…",3,[[["vec"],["usize"],["eta"]],[["partitionbound"],["str"],["result",["partitionbound","str"]]]]],[11,"from_noisy_estimates","","Partition bounds from noisy estimates",3,[[["vec"],["usize"],["eta"],["threadrandgen"],["option",["usize"]],["copy"],["partitionboundoptions"]],[["partitionbound"],["str"],["result",["partitionbound","str"]]]]],[0,"weights","b2dp::utilities","Compute weight tables via dynamic programming for the…",null,null],[3,"Key","b2dp::utilities::weights","",null,null],[12,"q","","Candidate value",4,null],[12,"i","","Index",4,null],[3,"WeightTable","","Weight table for sampling from the integer partition…",null,null],[12,"weights","","The weight table",5,null],[12,"arithmetic_config","","Arithmetic config",5,null],[12,"coefficients","","Coefficients table for bias computation",5,null],[12,"total_weights","","Total weights table for bias computation",5,null],[12,"probabilities","","Probabilities table for bias computation",5,null],[11,"from_bounds","","Computes the weight table for the given bounds, enforcing…",5,[[["eta"],["partitionbound"],["vec"]],[["str"],["result",["weighttable","str"]],["weighttable"]]]],[11,"inexact_from_bounds","","DO NOT USE except for bug demonstration! Computes the…",5,[[["eta"],["partitionbound"],["vec"]],[["str"],["result",["weighttable","str"]],["weighttable"]]]],[11,"get_bias","","Computes the bias of sampling from the weight table.…",5,[[["self"],["partitionbound"],["vec"]],[["vec",["f64"]],["str"],["result",["vec","str"]]]]],[0,"mechanisms","b2dp","Base-2 Differential Privacy Mechanisms",null,null],[0,"exponential","b2dp::mechanisms","Implements the base-2 exponential mechanism. ",null,null],[3,"ExponentialOptions","b2dp::mechanisms::exponential","The exponential mechanism optional parameters. ",null,null],[12,"min_retries","","The minimum number of retries for timing channel…",6,null],[12,"optimized_sample","","Whether to optimize sampling default: `false` Optimized…",6,null],[12,"empirical_precision","","Whether to use empirical precision default: `false`…",6,null],[5,"exponential_mechanism","","Implements the base-2 exponential mechanism. Utility…",null,[[["i64"],["u32"],["threadrandgen"],["exponentialoptions"],["copy"],["eta"],["fn"],["vec"]],[["t"],["result",["str"]],["str"]]]],[0,"laplace","b2dp::mechanisms","Implements an approximation of the Laplace mechanims via…",null,null],[5,"clamped_laplace_mechanism","b2dp::mechanisms::laplace","Implements clamped Laplace mechanism via the base-2…",null,[[["eta"],["threadrandgen"],["exponentialoptions"],["f64"],["copy"]],[["result",["f64","str"]],["f64"],["str"]]]],[0,"naive","b2dp::mechanisms","Naive Mechanism Implementations for comparison",null,null],[5,"naive_exponential_mechanism","b2dp::mechanisms::naive","",null,[[["vec"],["threadrandgen"],["f64"],["copy"]],[["t"],["result",["str"]],["str"]]]],[0,"sparsevector","b2dp::mechanisms","Implements the base-2 sparse vector mechanism.",null,null],[5,"sparse_vector","b2dp::mechanisms::sparsevector","The sparse vector mechanism Sensitivity (Delta) is assumed…",null,[[["bool"],["eta"],["usize"],["threadrandgen"],["vec"],["f64"],["copy"]],[["result",["vec","str"]],["vec",["bool"]],["str"]]]],[0,"integerpartition","b2dp::mechanisms","Implements the integer partition exponential mechanism…",null,null],[3,"IntegerPartitionOptions","b2dp::mechanisms::integerpartition","",null,null],[12,"optimize","","Whether to optimize sampling, range computation, etc.…",7,null],[5,"integer_partition_mechanism_with_weights","","The integer partition mechanism given a pre-computed…",null,[[["weighttable"],["partitionbound"],["integerpartitionoptions"]],[["str"],["vec",["i64"]],["result",["vec","str"]]]]],[5,"integer_partition_mechanism_with_bounds","","The integer partition mechanism given bounds",null,[[["eta"],["vec"],["partitionbound"],["integerpartitionoptions"]],[["str"],["vec",["i64"]],["result",["vec","str"]]]]],[11,"from","b2dp::utilities::params","",0,[[["t"]],["t"]]],[11,"into","","",0,[[],["u"]]],[11,"to_owned","","",0,[[["self"]],["t"]]],[11,"clone_into","","",0,[[["self"],["t"]]]],[11,"try_from","","",0,[[["u"]],["result"]]],[11,"try_into","","",0,[[],["result"]]],[11,"borrow","","",0,[[["self"]],["t"]]],[11,"borrow_mut","","",0,[[["self"]],["t"]]],[11,"type_id","","",0,[[["self"]],["typeid"]]],[11,"checked_as","","",0,[[],["option"]]],[11,"az","","",0,[[],["dst"]]],[11,"saturating_as","","",0,[[],["dst"]]],[11,"wrapping_as","","",0,[[],["dst"]]],[11,"overflowing_as","","",0,[[]]],[11,"static_as","","",0,[[],["option"]]],[11,"from","b2dp::utilities::exactarithmetic","",1,[[["t"]],["t"]]],[11,"into","","",1,[[],["u"]]],[11,"try_from","","",1,[[["u"]],["result"]]],[11,"try_into","","",1,[[],["result"]]],[11,"borrow","","",1,[[["self"]],["t"]]],[11,"borrow_mut","","",1,[[["self"]],["t"]]],[11,"type_id","","",1,[[["self"]],["typeid"]]],[11,"checked_as","","",1,[[],["option"]]],[11,"az","","",1,[[],["dst"]]],[11,"saturating_as","","",1,[[],["dst"]]],[11,"wrapping_as","","",1,[[],["dst"]]],[11,"overflowing_as","","",1,[[]]],[11,"static_as","","",1,[[],["option"]]],[11,"from","b2dp::utilities::randomness","",8,[[["t"]],["t"]]],[11,"into","","",8,[[],["u"]]],[11,"to_owned","","",8,[[["self"]],["t"]]],[11,"clone_into","","",8,[[["self"],["t"]]]],[11,"try_from","","",8,[[["u"]],["result"]]],[11,"try_into","","",8,[[],["result"]]],[11,"borrow","","",8,[[["self"]],["t"]]],[11,"borrow_mut","","",8,[[["self"]],["t"]]],[11,"type_id","","",8,[[["self"]],["typeid"]]],[11,"checked_as","","",8,[[],["option"]]],[11,"az","","",8,[[],["dst"]]],[11,"saturating_as","","",8,[[],["dst"]]],[11,"wrapping_as","","",8,[[],["dst"]]],[11,"overflowing_as","","",8,[[]]],[11,"static_as","","",8,[[],["option"]]],[11,"from","b2dp::utilities::bounds","",2,[[["t"]],["t"]]],[11,"into","","",2,[[],["u"]]],[11,"try_from","","",2,[[["u"]],["result"]]],[11,"try_into","","",2,[[],["result"]]],[11,"borrow","","",2,[[["self"]],["t"]]],[11,"borrow_mut","","",2,[[["self"]],["t"]]],[11,"type_id","","",2,[[["self"]],["typeid"]]],[11,"checked_as","","",2,[[],["option"]]],[11,"az","","",2,[[],["dst"]]],[11,"saturating_as","","",2,[[],["dst"]]],[11,"wrapping_as","","",2,[[],["dst"]]],[11,"overflowing_as","","",2,[[]]],[11,"static_as","","",2,[[],["option"]]],[11,"from","","",3,[[["t"]],["t"]]],[11,"into","","",3,[[],["u"]]],[11,"try_from","","",3,[[["u"]],["result"]]],[11,"try_into","","",3,[[],["result"]]],[11,"borrow","","",3,[[["self"]],["t"]]],[11,"borrow_mut","","",3,[[["self"]],["t"]]],[11,"type_id","","",3,[[["self"]],["typeid"]]],[11,"checked_as","","",3,[[],["option"]]],[11,"az","","",3,[[],["dst"]]],[11,"saturating_as","","",3,[[],["dst"]]],[11,"wrapping_as","","",3,[[],["dst"]]],[11,"overflowing_as","","",3,[[]]],[11,"static_as","","",3,[[],["option"]]],[11,"from","b2dp::utilities::weights","",4,[[["t"]],["t"]]],[11,"into","","",4,[[],["u"]]],[11,"to_owned","","",4,[[["self"]],["t"]]],[11,"clone_into","","",4,[[["self"],["t"]]]],[11,"try_from","","",4,[[["u"]],["result"]]],[11,"try_into","","",4,[[],["result"]]],[11,"borrow","","",4,[[["self"]],["t"]]],[11,"borrow_mut","","",4,[[["self"]],["t"]]],[11,"type_id","","",4,[[["self"]],["typeid"]]],[11,"checked_as","","",4,[[],["option"]]],[11,"az","","",4,[[],["dst"]]],[11,"saturating_as","","",4,[[],["dst"]]],[11,"wrapping_as","","",4,[[],["dst"]]],[11,"overflowing_as","","",4,[[]]],[11,"static_as","","",4,[[],["option"]]],[11,"from","","",5,[[["t"]],["t"]]],[11,"into","","",5,[[],["u"]]],[11,"try_from","","",5,[[["u"]],["result"]]],[11,"try_into","","",5,[[],["result"]]],[11,"borrow","","",5,[[["self"]],["t"]]],[11,"borrow_mut","","",5,[[["self"]],["t"]]],[11,"type_id","","",5,[[["self"]],["typeid"]]],[11,"checked_as","","",5,[[],["option"]]],[11,"az","","",5,[[],["dst"]]],[11,"saturating_as","","",5,[[],["dst"]]],[11,"wrapping_as","","",5,[[],["dst"]]],[11,"overflowing_as","","",5,[[]]],[11,"static_as","","",5,[[],["option"]]],[11,"from","b2dp::mechanisms::exponential","",6,[[["t"]],["t"]]],[11,"into","","",6,[[],["u"]]],[11,"to_owned","","",6,[[["self"]],["t"]]],[11,"clone_into","","",6,[[["self"],["t"]]]],[11,"try_from","","",6,[[["u"]],["result"]]],[11,"try_into","","",6,[[],["result"]]],[11,"borrow","","",6,[[["self"]],["t"]]],[11,"borrow_mut","","",6,[[["self"]],["t"]]],[11,"type_id","","",6,[[["self"]],["typeid"]]],[11,"checked_as","","",6,[[],["option"]]],[11,"az","","",6,[[],["dst"]]],[11,"saturating_as","","",6,[[],["dst"]]],[11,"wrapping_as","","",6,[[],["dst"]]],[11,"overflowing_as","","",6,[[]]],[11,"static_as","","",6,[[],["option"]]],[11,"from","b2dp::mechanisms::integerpartition","",7,[[["t"]],["t"]]],[11,"into","","",7,[[],["u"]]],[11,"to_owned","","",7,[[["self"]],["t"]]],[11,"clone_into","","",7,[[["self"],["t"]]]],[11,"try_from","","",7,[[["u"]],["result"]]],[11,"try_into","","",7,[[],["result"]]],[11,"borrow","","",7,[[["self"]],["t"]]],[11,"borrow_mut","","",7,[[["self"]],["t"]]],[11,"type_id","","",7,[[["self"]],["typeid"]]],[11,"checked_as","","",7,[[],["option"]]],[11,"az","","",7,[[],["dst"]]],[11,"saturating_as","","",7,[[],["dst"]]],[11,"wrapping_as","","",7,[[],["dst"]]],[11,"overflowing_as","","",7,[[]]],[11,"static_as","","",7,[[],["option"]]],[11,"clone","b2dp::utilities::params","",0,[[["self"]],["eta"]]],[11,"clone","b2dp::utilities::randomness","",8,[[["self"]],["generatoropenssl"]]],[11,"clone","b2dp::utilities::weights","",4,[[["self"]],["key"]]],[11,"clone","b2dp::mechanisms::exponential","",6,[[["self"]],["exponentialoptions"]]],[11,"clone","b2dp::mechanisms::integerpartition","",7,[[["self"]],["integerpartitionoptions"]]],[11,"default","b2dp::utilities::bounds","Default options for generating PartitionBounds…",2,[[],["partitionboundoptions"]]],[11,"default","b2dp::mechanisms::exponential","Default options for the exponential mechanism `min_retries…",6,[[],["exponentialoptions"]]],[11,"default","b2dp::mechanisms::integerpartition","Default options for the exponential mechanism `min_retries…",7,[[],["integerpartitionoptions"]]],[11,"eq","b2dp::utilities::params","",0,[[["self"]],["bool"]]],[11,"eq","b2dp::utilities::weights","",4,[[["self"],["key"]],["bool"]]],[11,"ne","","",4,[[["self"],["key"]],["bool"]]],[11,"fmt","b2dp::utilities::params","",0,[[["formatter"],["self"]],["result"]]],[11,"fmt","b2dp::utilities::exactarithmetic","",1,[[["formatter"],["self"]],["result"]]],[11,"fmt","b2dp::utilities::randomness","",8,[[["formatter"],["self"]],["result"]]],[11,"fmt","b2dp::utilities::bounds","",2,[[["formatter"],["self"]],["result"]]],[11,"fmt","","",3,[[["formatter"],["self"]],["result"]]],[11,"fmt","b2dp::utilities::weights","",4,[[["formatter"],["self"]],["result"]]],[11,"fmt","","",5,[[["formatter"],["self"]],["result"]]],[11,"fmt","b2dp::mechanisms::exponential","",6,[[["formatter"],["self"]],["result"]]],[11,"fmt","b2dp::mechanisms::integerpartition","",7,[[["formatter"],["self"]],["result"]]],[11,"hash","b2dp::utilities::weights","",4,[[["self"],["__h"]]]],[11,"gen","b2dp::utilities::randomness","",8,[[["self"]],["u32"]]]],"p":[[3,"Eta"],[3,"ArithmeticConfig"],[3,"PartitionBoundOptions"],[3,"PartitionBound"],[3,"Key"],[3,"WeightTable"],[3,"ExponentialOptions"],[3,"IntegerPartitionOptions"],[3,"GeneratorOpenSSL"]]};
addSearchOptions(searchIndex);initSearch(searchIndex);