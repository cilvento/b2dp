//! Implements the base-2 sparse vector mechanism.

use rug::{Float, rand::ThreadRandGen};
use crate::utilities::exactarithmetic::{ArithmeticConfig};
use crate::utilities::params::Eta;
use crate::utilities::discretesampling::{sample_within_bounds, noisy_threshold, adjust_eta, is_multiple_of};
use crate::errors::*;

/// The sparse vector mechanism
/// Sensitivity (Delta) is assumed to be 1. Privacy parameters must be scaled 
/// appropriately if this is not the case. 
/// 
/// ## Arguments:
///   * `eta1`: the privacy budget to spend on sampling the threshold
///   * `eta2`: the privacy budget to spend on the queries. Sensitivity of 
///         queries assumed to be 1. `eta2` must be scaled appropriately for c.
///   * `c`: the number of allowed positive queries. 
///   * `gamma`: the granularity (must be a reciprocal of a positive integer)
///   * `Q`: a set of query values which must be integer multiples of `gamma`, otherwise
///          they are rounded, which must be accounted for in privacy budget calculations.
///   * `Qmin`: the minimum value any query may be.
///   * `Qmax`: the maximum value any query may be. 
///   * `w`: the width to use for threshold comparison (higher width increases 
///         accuracy at the cost of efficiency, but only up to a point)
///   * `rng`: randomness source
///   * `optimize`: whether to optimize sampling (exacerbates timing channels)
/// 
/// ## Privacy budget usage
/// Uses `eta1*Delta + 2*eta2*Delta*c`. As such, `eta` and `eta2` should be 
/// scaled appropriately for the desired usage. 
/// 
/// ## Returns
/// A vector of `bool` with at most `c` `true` entries indicating which queries 
/// exceeded the noisy threshold. 
/// 
/// ## Exact Arithmetic
/// This method calls `enter_exact_scope` which clears the `mpfr::flags`. 
/// 
/// ## Timing Channels
/// This method uses [normalized_sample](../../utilities/exactarithmetic/fn.normalized_sample.html#known-timing-channels)
/// which has known timing channels. Care should be taken given that 
/// the method allows for an arbitrary number of queries which may amplify
/// timing channels. 
pub fn sparse_vector<R: ThreadRandGen+Copy>(eta1: Eta, eta2: Eta, 
                                            c: usize, query_values: &Vec<f64>, 
                                            gamma: f64, 
                                            query_min: f64, query_max: f64, 
                                            w: f64, rng: R,
                                            optimize: bool)
    -> Result<Vec<bool>>
{
    // Initialize output vector
    let mut outputs: Vec<bool> = Vec::new();
    let mut count = 0;
    
    // Get an arithmetic_config to keep track of required precision
    let mut arithmetic_config = ArithmeticConfig::basic()?;
    // enter exact scope
    arithmetic_config.enter_exact_scope()?;
    // Check gamma and adjustments to ensure sufficient precision
    // `adjust_eta` also checks for validity of gamma and that
    // g_inv is integer. 
    let mut g = arithmetic_config.get_float(gamma);
    let mut g_inv = arithmetic_config.get_float(1.0/gamma);
    let _eta1_prime = adjust_eta(eta1, &g_inv, & mut arithmetic_config)?; 
    let _eta2_prime = adjust_eta(eta2, &g_inv, & mut arithmetic_config)?;

    // If inexact arithmetic performed, increase precision and try again
    while arithmetic_config.exit_exact_scope().is_err() {
        arithmetic_config.increase_precision(16)?; // increase precision
        arithmetic_config.inexact_arithmetic = false; // reset the inexact arithmetic flag on the configuration
        arithmetic_config.enter_exact_scope()?; // re-enter the exact scope, and try again
        g = arithmetic_config.get_float(gamma);
        g_inv = arithmetic_config.get_float(1.0/gamma);
        let _eta1_prime = adjust_eta(eta1, &g_inv, & mut arithmetic_config)?;
        let _eta2_prime = adjust_eta(eta2, &g_inv, & mut arithmetic_config)?;
    }
    // re-enter the exact scope after appropriate precision determined
    arithmetic_config.enter_exact_scope()?;

    // Check that w, Qmax, Qmin, etc are multiples of gamma
    if !is_multiple_of(&arithmetic_config.get_float(w),&g) 
        { return Err("w is not an integer multiple of gamma.".into()); }
    if !is_multiple_of(&arithmetic_config.get_float(query_min),&g) 
        { return Err("query_min is not an integer multiple of gamma.".into()); }
    if !is_multiple_of(&arithmetic_config.get_float(query_max),&g) 
        { return Err("query_max is not an integer multiple of gamma.".into()); }
    

    // Sample Rho
    let rho = sample_within_bounds(eta1,
                                    &arithmetic_config.get_float(&g),
                                    &arithmetic_config.get_float(query_min - w),
                                    &arithmetic_config.get_float(query_max + w),
                                    &mut arithmetic_config,
                                    rng,
                                    optimize)?;
    // Iterate through queries 
    for i in 0..query_values.len() {
        // Clamp Q[i] if needed
        let mut q: Float;
        if query_values[i] > query_max { q = arithmetic_config.get_float(query_max); }
        else if query_values[i] < query_min { q = arithmetic_config.get_float(query_min); }
        else { q = arithmetic_config.get_float(query_values[i]); }
        // Check that q is a multiple of gamma
        if !is_multiple_of(&q,&g) { 
            // Round
            let mut m = arithmetic_config.get_float(&q/&g);
            m.round_mut(); 
            q = arithmetic_config.get_float(&g*&m);
        }
    

        // Compute rho hat
        let rho_hat: Float;
        let rho_max = arithmetic_config.get_float(&q + w);
        let rho_min = arithmetic_config.get_float(&q - w);
        if rho > rho_max { rho_hat = rho_max; }
        else if rho < rho_min {rho_hat = rho_min; } 
        else { rho_hat = arithmetic_config.get_float(&rho); }

        // Run noisy threshold
        let g_inv = arithmetic_config.get_float(1.0/gamma);
        let a = noisy_threshold(eta2,& mut arithmetic_config,&g_inv,&rho_hat,rng,optimize)?;
        if a.is_infinite() && a.is_sign_positive() { 
            outputs.push(true); 
            count += 1;
            if count >= c { return Ok(outputs); } // if we have already encountered c positives, stop
        }
        else if a.is_infinite() && a.is_sign_negative() { outputs.push(false); }
    }
    return Ok(outputs);
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::utilities::randomness::GeneratorOpenSSL;


    #[test]
    fn test_sparse_vector(){
        let eta1 = Eta::new(1,1,2).unwrap();
        let eta2 = Eta::new(1,1,2).unwrap();
        let c = 2;
        let queries = vec![1.0,2.0,3.0,4.0,5.0,1.0];
        let gamma = 0.5;
        let q_min = 0.0;
        let q_max = 6.0;
        let w = 5.0;
        let rng = GeneratorOpenSSL {};
        let optimize = false;
        let outputs = sparse_vector(eta1, eta2, c, &queries, gamma, q_min, q_max, w, rng, optimize);
        assert!(outputs.is_ok());
        println!("{:?}", outputs);
    }
}