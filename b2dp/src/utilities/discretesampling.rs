//! Implements discrete sampling methods in base-2
//! including [`noisy_threshold`](./fn.noisy_threshold.html) and 
//! [`sample_within_bounds`](./fn.sample_within_bounds.html).

use rug::{ops::Pow, rand::ThreadRandGen, Float, float::Special};
use super::exactarithmetic::{ArithmeticConfig, normalized_sample};
use super::params::Eta;
use crate::errors::*;
/// Returns whether a is an integer multiple of b
pub fn is_multiple_of(a: &Float, b: &Float) -> bool {
    if a.is_infinite() || b.is_infinite() { return false; }
    let c = Float::with_val(a.prec(),a).remainder(b);
    if c == 0 {
        let d = Float::with_val(a.prec(), a/b);
        return d.is_integer();
    }
    return false;
}

/// Samples from the Discrete Laplace mechanism at granularity
/// `gamma` within the provided bounds `wmin` and `wmax` where
/// the values `wmin` and `wmax` are sampled with probability
/// equal to the sum of the probabilities of all elements
/// greater than or less than the bounds, respectively. 
/// ## Arguments 
///  * `eta`: privacy parameter
///  * `gamma`: granularity parameter, must be reciprocal of an integer. 
///  * `wmin`: the minimum bound, must be an integer multiple of `gamma`.
///  * `wmax`: the maximum bound, must be an integer multiple of `gamma`.
///  * `arithmeticconfig`: ArithmeticConfig with sufficient precision for
///      sampling. 
///  * `rng`: randomness source
///  * `optimize`: whether to optimize sampling. 
/// 
/// ## Returns
/// Returns an integer multiple of `gamma` within the provided bounds
/// sampled according to the discrete laplace mechanism or an error if
/// `eta` cannot be adjusted by `gamma`.
/// 
/// ## Privacy Budget Usage
/// Uses `eta` privacy budget.  
/// 
/// ## Exact Arithmetic
/// Does not enforce exact arithmetic, this is the caller's responsibility. 
/// 
/// ## Timing Channels
/// Uses [`normalized_sample`](../exactarithmetic/fn.normalized_sample.html#known-timing-channels) 
/// which has known timing channels. 
/// 
/// ## Example Usage 
/// ```
/// # use b2dp::{Eta,GeneratorOpenSSL,utilities::exactarithmetic::ArithmeticConfig, sample_within_bounds};
/// # use rug::Float;
/// # use b2dp::errors::*;
/// # fn main() -> Result<()> {
/// let eta = Eta::new(1,1,2)?; // construct eta that can be adjusted for the desired value of gamma.
/// let mut arithmeticconfig = ArithmeticConfig::basic()?;
/// let rng = GeneratorOpenSSL {};
/// let gamma = arithmeticconfig.get_float(0.5);
/// let wmin = arithmeticconfig.get_float(-5);
/// let wmax = arithmeticconfig.get_float(5);
/// arithmeticconfig.enter_exact_scope()?;
/// let s = sample_within_bounds(eta, &gamma, &wmin, &wmax, & mut arithmeticconfig, rng,false)?;
/// let b = arithmeticconfig.exit_exact_scope();
/// assert!(b.is_ok()); // Must check that no exact arithmetic was performed. 
/// # Ok(())
/// # }
/// ```
pub fn sample_within_bounds<R: ThreadRandGen+Copy>(eta: Eta, gamma: &Float, 
                                                wmin: &Float, wmax: &Float, 
                                                arithmeticconfig: & mut ArithmeticConfig,
                                                rng: R, optimize: bool)
    -> Result<Float>
{
    // Check inputs
    if wmax <= wmin { return Err("`wmin` must be strictly less than `wmax`.".into()); }
    if !is_multiple_of(&wmin, &gamma) { return Err("`wmin` is not integer multiple of `gamma`.".into()); }
    if !is_multiple_of(&wmax, &gamma) { return Err("`wmax` is not integer multiple of `gamma`.".into()); }
    let t_min = arithmeticconfig.get_float(wmin/gamma);
    let t_max = arithmeticconfig.get_float(wmax/gamma);

    // Adjust eta
    let gamma_inv = arithmeticconfig.get_float(1/gamma);
    let eta_prime = adjust_eta(eta, &gamma_inv, arithmeticconfig)?;
    let base = eta_prime.get_base(arithmeticconfig.precision)?;
    let plus_infty = arithmeticconfig.get_float(Special::Infinity);
    let neg_infty = arithmeticconfig.get_float(Special::NegInfinity);

    // Get the weights for each region
    let p_l = get_sum(&base, arithmeticconfig, 
                        &neg_infty,
                        &t_min)?;
    let p_u = get_sum(&base, arithmeticconfig, 
                        &t_max,
                        &plus_infty)?;
    let p_t = get_sum(&base, arithmeticconfig, 
                        &neg_infty,
                        &plus_infty)?;
    let p_m = arithmeticconfig.get_float(p_t - &p_u - &p_l);
    
    let region_weights: Vec<Float> = vec![p_l,p_u,p_m];

    // Sample which region rho lies in
    let r = normalized_sample(&region_weights,arithmeticconfig,rng,optimize)?;

    match r { 
        0 => return Ok(arithmeticconfig.get_float(Special::NegInfinity)),  // lower region
        1 => return Ok(arithmeticconfig.get_float(Special::Infinity)),     // upper region
        _=> (), // Must sample from the middle region 
    }
    // Sample from the middle region
    // Construct weights and outcome space
    let mut outcomes: Vec<Float> = Vec::new();
    let mut weights: Vec<Float> = Vec::new();
    let mut k = arithmeticconfig.get_float(wmin/gamma);
    k = k+1;// 
    let mut o = arithmeticconfig.get_float(&k*gamma);

    while o < *wmax {
        // record the outcome and the weight
        outcomes.push(o);
        let next_k = arithmeticconfig.get_float(&k+1);
        let w = arithmeticconfig.get_float((&base).pow(k.abs()));
        weights.push(w);
        // increment k and update the outcome
        k = next_k;
        o = arithmeticconfig.get_float(&k*gamma);
    } 

    // Sample
    let s = normalized_sample(&weights,arithmeticconfig,rng,optimize)?;
    return Ok(arithmeticconfig.get_float(&outcomes[s]));
}


/// Adjust the given eta to account for granularity 1/gamma_inv
/// 
/// ## Arguments
/// * `eta`: privacy parameter
/// * `gamma_inv`: an integer-valued Float representing the inverse of `gamma`
/// * `arithmeticconfig`: an ArithmeticConfig with sufficient precision to adjust
///    `eta` if it can be adjusted. 
/// 
/// ## Returns
/// A new `Eta` adjusted such that `2^eta_prime = (2^eta)^gamma` or an error
/// if `eta` cannot be adjusted. 
/// 
/// ## Exact Arithmetic
/// Does not enforce exact arithmetic, this is the caller's responsibility.
pub fn adjust_eta(eta: Eta, gamma_inv: &Float, arithmeticconfig: & mut ArithmeticConfig) 
    -> Result<Eta>//,&'static str> 
{
    // Check that gamma is valid for the given eta
    if !gamma_inv.is_integer() { return Err("`gamma_inv` must be an integer.".into()); }
    let gamma = arithmeticconfig.get_float(1.0/gamma_inv);
    // Check if eta.z is divisible by gamma_inv.
    let mut z_prime = eta.z; 
    let mut x_prime = eta.x;
    let mut y_prime = eta.y;
    if arithmeticconfig.get_float(eta.z*&gamma).is_integer()
    {
        let rootz = arithmeticconfig.get_float(eta.z*&gamma);
        
        z_prime = rootz.to_integer().unwrap().to_u32().unwrap(); 
        // Leave x and y as is
    }
    else { 
        // Leave z as is
        // Check if x and y meet the critera
        let fx = arithmeticconfig.get_float(eta.x);
        let fy = arithmeticconfig.get_float(eta.y);
        
        if !is_multiple_of(&fy, &gamma) {return Err("Unable to adjust for gamma (y).".into());}
        let rooty = arithmeticconfig.get_float(fy*&gamma);
        
        let rootx = arithmeticconfig.get_float(fx.pow(&gamma));
        if !rootx.is_integer() {return Err("Unable to adjust for gamma (x).".into());}
        
        x_prime = rootx.to_integer().unwrap().to_u32().unwrap(); // TODO: more elegant error handling
        y_prime = rooty.to_integer().unwrap().to_u32().unwrap();
    }
    let eta_prime = Eta::new(x_prime,y_prime,z_prime)?;
    return Ok(eta_prime);

}


/// Determines whether the discretized Laplace exceeds the given threshold
/// ## Arguments
///   * `eta`: the privacy parameter
///   * `arithmeticconfig`:
///   * `gamma`
///   * `threshold`: 
///   * `rng`: randomness source
///   * `optimize`: whether to optimize sampling, exacerbates timing channels
/// 
/// ## Returns
/// Returns a `Float` with value `Special::Infinity` if greater than or equal to the threshold,
/// otherwise returns with value `Special::NegInfinity`. Returns an error if `eta` cannot
/// be appropriately adjusted or sum computation fails. 
/// 
/// ## Exact Arithmetic
/// Does not explicitly enforce exact arithmetic, this is the caller's responsibility.
/// 
/// ## Privacy Budget
/// Uses `eta` privacy budget
/// 
/// ## Timing Channels
/// Uses [`normalized_sample`](../exactarithmetic/fn.normalized_sample.html#known-timing-channels) 
/// which has known timing channels. 
/// 
/// ## Example Usage
/// ```
/// # use b2dp::{Eta,GeneratorOpenSSL,utilities::exactarithmetic::ArithmeticConfig, noisy_threshold};
/// # use rug::Float;
/// # use b2dp::errors::*;
/// # fn main() -> Result<()> {
/// let eta = Eta::new(1,1,2)?; // construct eta that can be adjusted for the desired value of gamma.
/// let mut arithmeticconfig = ArithmeticConfig::basic()?;
/// let rng = GeneratorOpenSSL {};
/// let gamma_inv = arithmeticconfig.get_float(2);
/// let threshold = arithmeticconfig.get_float(0);
/// arithmeticconfig.enter_exact_scope()?; 
/// let s = noisy_threshold(eta, & mut arithmeticconfig, &gamma_inv, &threshold, rng, false)?;
/// assert!(!s.is_finite()); // returns plus or minus infinity
/// if s.is_sign_positive() { /* Greater than the threshold */ ;}
/// else { /* Less than the threshold. */ ;}
/// let b = arithmeticconfig.exit_exact_scope();
/// assert!(b.is_ok()); // Must check that no exact arithmetic was performed. 
/// # Ok(())
/// # }
/// ```
pub fn noisy_threshold<R: ThreadRandGen>(eta: Eta, arithmeticconfig: & mut ArithmeticConfig, gamma_inv: &Float, threshold: &Float,
                        rng: R, optimize: bool) 
    -> Result<Float>
    { 
        // plus and minus infinity
        let plus_infty = arithmeticconfig.get_float(Special::Infinity);
        let neg_infty = arithmeticconfig.get_float(Special::NegInfinity);

        // Adjust eta to take gamma into account.
        let eta_prime = adjust_eta(eta, &gamma_inv, arithmeticconfig)?;
        // get the base for the adjusted eta
        let base = eta_prime.get_base(arithmeticconfig.precision)?;
        
        // Check that gamma is valid (integer reciprocal)
        if !gamma_inv.is_integer() { return Err("`gamma_inv` must be an integer.".into()); }
        let gamma = arithmeticconfig.get_float(1/gamma_inv);

        // Check that threshold is integer multiple of gamma
        if !is_multiple_of(&threshold, &gamma)  { return Err("`threshold` must be integer multiple of `gamma`.".into()); }
        let t = arithmeticconfig.get_float(threshold*gamma_inv); 

        let p_top = get_sum(&base,
                            arithmeticconfig, 
                            &t, 
                            &plus_infty)?;
        let p_total = get_sum(&base,
                            arithmeticconfig, 
                            &neg_infty, 
                            &plus_infty)?;
        let p_bot = arithmeticconfig.get_float(p_total - &p_top);
        let weights: Vec<Float> = vec![p_top,p_bot];
        let s = normalized_sample(&weights, arithmeticconfig,rng,optimize)?;
        
        if s == 0 {
            return Ok(arithmeticconfig.get_float(Special::Infinity));
        }
        else {
            return Ok(arithmeticconfig.get_float(Special::NegInfinity));
        }
    }


/// Returns the sum: 
/// 0.5*(1-base)^{-1}*\sum_{k=start}^{end}base^{|k|}
/// 
/// ## Arguments:
///   * `base`: a `Float` indicating the base for the sum
///   * `arithmeticconfig`: an ArithmeticConfig  
///   * `start`: an integer-valued or infinite `Float` indicating the starting point of the sum
///   * `end`: an integer-valued or infinite `Float` indicating the ending point of the sum
/// 
/// ## Returns
/// Either the sum or an error if parameters are mis-specified. This method does not explicitly check
/// for inexact arithmetic, it is the caller's responsibility to do so. 
/// 
/// ## Exact Arithmetic
/// This method does not enforce exact arithmetic, this is the caller's responsibility. 
/// 
/// ## Timing channels
/// The recursive calls to `get_sum` introduce a timing channel distinguishing whether
/// `start` and `end` cross zero. 
pub fn get_sum(base: &Float, arithmeticconfig: &ArithmeticConfig, start: &Float, end: &Float) 
    -> Result<Float>
{
    // Check ordering
    if start >= end { return Err("`start` must be strictly less than `end`.".into()); }
    // Check integrality
    if (!start.is_integer() && start.is_finite()) || // infinite values are not considered integers
        (!end.is_integer() && end.is_finite()) { 
        return Err("`start` and `end` must be integer values or infinite.".into());
    }

    // Check base magnitude
    if *base >= 1 { return Err("`base` must be less than 1.".into()); }

    // Check for negative infinity case
    if start.is_infinite() && start.is_sign_negative() {

        if end.is_infinite() && end.is_sign_positive() {
            // Full infinite sum 
            // = 0.5*(1+base)
            let base_plus_one = arithmeticconfig.get_float(1+base); 
            return Ok(arithmeticconfig.get_float(0.5 * base_plus_one));
        }

        // Half-open infinite sum
        let m = arithmeticconfig.get_float(end);
        let abs_end = m.abs();
        if *end < 0 {
            // Same sign
            // = 0.5 - 0.5*(1-base^(|end|)) 
            let p = arithmeticconfig.get_float(base).pow(abs_end); // base^(|end|)
            let s = arithmeticconfig.get_float(0.5 - 0.5*(1-p));
            return Ok(s);
        }
        else {
            // Different sign
            // = 0.5 + 0.5*base - 0.5*base^(|end| + 1)
            let end_plus_one = arithmeticconfig.get_float(abs_end + 1);
            let p = arithmeticconfig.get_float(base).pow(end_plus_one); // base^(|end|+1)
            let half_base = arithmeticconfig.get_float(0.5*base);
            let s = arithmeticconfig.get_float(0.5 + half_base - 0.5*p);
            return Ok(s);
        }

    }
    // Half-open positive infinite sum
    else if end.is_infinite() && end.is_sign_positive() {
        let m = arithmeticconfig.get_float(start);
        let abs_start = m.abs();
    
        // Half-open infinite sum
        if *start > 0 {
            // Same sign
            // = 0.5 - 0.5*(1-base^(|start|))
            let p = arithmeticconfig.get_float(base).pow(abs_start); // base^(|start|)
           
            
            let s = arithmeticconfig.get_float(0.5 - 0.5*(1-p));
            return Ok(s);
        }
        else {
            // Different sign
            // = 0.5 + 0.5*base - 0.5*base^(|start| + 1)
            let start_plus_one = arithmeticconfig.get_float(abs_start + 1);
           
            let p = arithmeticconfig.get_float(base).pow(start_plus_one); // base^(|start|+1)
           
            let half_base = arithmeticconfig.get_float(0.5*base);
            let s = arithmeticconfig.get_float(0.5 + half_base - 0.5*p);
            return Ok(s);
        }
    }
    
    
    // Otherwise, finite sum, recurse 
    // = get_sum(-infinity,infinity) - get_sum(-infinity, start - 1) - get_sum(end + 1, infinity)
    let plus_infty = arithmeticconfig.get_float(Special::Infinity);
    let neg_infty = arithmeticconfig.get_float(Special::NegInfinity);
    
    let total_sum = get_sum(base,
                            arithmeticconfig,
                            &neg_infty,
                            &plus_infty)?;
    let neg_sum = get_sum(base,
                            arithmeticconfig, 
                            &neg_infty, 
                            &arithmeticconfig.get_float(start-1))?;
    let pos_sum = get_sum(base,
                            arithmeticconfig,
                            &arithmeticconfig.get_float(end+1), 
                            &plus_infty)?;
    
    let s = arithmeticconfig.get_float(total_sum - neg_sum - pos_sum);
    Ok(s)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::utilities::randomness::GeneratorOpenSSL;


    #[test]
    fn test_sample_within_bounds(){
        let eta = Eta::new(1,1,2).unwrap();
        let rng = GeneratorOpenSSL {};
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let gamma = arithmeticconfig.get_float(0.5);
        let gamma_inv = arithmeticconfig.get_float(1/&gamma);
        let eta_prime = adjust_eta(eta, &gamma_inv, &mut arithmeticconfig);
        assert!(eta_prime.is_ok());
        let wmin = arithmeticconfig.get_float(-5);
        let wmax = arithmeticconfig.get_float(5);
        let _a = arithmeticconfig.enter_exact_scope();
        let s = sample_within_bounds(eta, &gamma, &wmin, &wmax, & mut arithmeticconfig, rng,false);
        
        assert!(s.is_ok());
        let b = arithmeticconfig.exit_exact_scope();
        assert!(b.is_ok());
    }
    #[test]
    fn test_sample_within_bounds_inf_likely(){
        let eta = Eta::new(1,1,2).unwrap();
        let rng = GeneratorOpenSSL {};
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let gamma = arithmeticconfig.get_float(0.5);
        let gamma_inv = arithmeticconfig.get_float(1/&gamma);
        let eta_prime = adjust_eta(eta, &gamma_inv, &mut arithmeticconfig);
        assert!(eta_prime.is_ok());
        let wmin = arithmeticconfig.get_float(-1);
        let wmax = arithmeticconfig.get_float(1);
        let _a = arithmeticconfig.enter_exact_scope();
        let s = sample_within_bounds(eta, &gamma, &wmin, &wmax, & mut arithmeticconfig, rng,false);
        
        assert!(s.is_ok());
        let b = arithmeticconfig.exit_exact_scope();
        assert!(b.is_ok());
    }
    #[test]
    fn test_remainder(){
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let gamma = arithmeticconfig.get_float(0.25);
        let t = arithmeticconfig.get_float(0.75);
        let _a = arithmeticconfig.enter_exact_scope();
        let r = t.remainder(&gamma);
        println!("{:?}", &r);
        let b = arithmeticconfig.exit_exact_scope();
        assert!(b.is_ok());
        
    }


    #[test]
    fn test_noisy_threshold(){

        // Simple zero threshold test
        let eta = Eta::new(1,1,2).unwrap();
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let rng = GeneratorOpenSSL {};
        let gamma_inv = arithmeticconfig.get_float(2);
        let threshold = arithmeticconfig.get_float(0);
        let _a = arithmeticconfig.enter_exact_scope();
        let s = noisy_threshold(eta, & mut arithmeticconfig, &gamma_inv, &threshold, rng, false).unwrap();
        assert!(!s.is_finite()); // should get plus or minus infinity
        let b = arithmeticconfig.exit_exact_scope();
        assert!(b.is_ok());

        // Check fail on threshold that is not a multiple of gamma
        let eta = Eta::new(1,1,2).unwrap();
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let rng = GeneratorOpenSSL {};
        let gamma_inv = arithmeticconfig.get_float(2);
        let threshold = arithmeticconfig.get_float(0.3);
        let _a = arithmeticconfig.enter_exact_scope();
        let s = noisy_threshold(eta, & mut arithmeticconfig, &gamma_inv, &threshold, rng, false);
        assert!(s.is_err());
        let b = arithmeticconfig.exit_exact_scope();
        assert!(b.is_ok());
        
        // Check fail on eta that cannot be adjusted
        let eta = Eta::new(1,1,1).unwrap();
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let rng = GeneratorOpenSSL {};
        let gamma_inv = arithmeticconfig.get_float(2);
        let threshold = arithmeticconfig.get_float(0.3);
        let _a = arithmeticconfig.enter_exact_scope();
        let s = noisy_threshold(eta, & mut arithmeticconfig, &gamma_inv, &threshold, rng, false);
        assert!(s.is_err());
        let _b = arithmeticconfig.exit_exact_scope(); 
    }

    // Test eta adjustment
    #[test]
    fn test_eta_adjustment(){
        // Test case passes by modifying z
        let eta = Eta::new(1,1,2).unwrap();
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let gamma_inv = arithmeticconfig.get_float(2);
        let eta_prime = adjust_eta(eta, &gamma_inv, & mut arithmeticconfig).unwrap();
        assert_eq!(eta_prime, Eta::new(1,1,1).unwrap());

        // Test case passes by modifying x and y
        let eta = Eta::new(1,2,1).unwrap();
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let gamma_inv = arithmeticconfig.get_float(2);
        let eta_prime = adjust_eta(eta, &gamma_inv, & mut arithmeticconfig).unwrap();
        assert_eq!(eta_prime, Eta::new(1,1,1).unwrap());

        // Cannot be adjusted
        let eta = Eta::new(1,1,1).unwrap();
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let gamma_inv = arithmeticconfig.get_float(2);
        let eta_prime = adjust_eta(eta, &gamma_inv, & mut arithmeticconfig);
        assert!(eta_prime.is_err());
    }

    // Test sum parameter errors
    #[test]
    fn test_sum_parameters(){
        let arithmeticconfig = ArithmeticConfig::basic().unwrap();

        // base that is too large
        let base = arithmeticconfig.get_float(1.0);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 0), &Float::with_val(32, 5));
        assert!(s.is_err());

        // start > end
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 5), &Float::with_val(32, 0));
        assert!(s.is_err());
        
        // start = end
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 5), &Float::with_val(32, 5));
        assert!(s.is_err());

        // infinite and equal but different precision
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity), &Float::with_val(16, Special::NegInfinity));
        assert!(s.is_err());

        // non-integer
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity), &Float::with_val(16, 1.75));
        assert!(s.is_err());

    }

    // Test Sum computation
    #[test]
    fn test_sums() {
        let arithmeticconfig = ArithmeticConfig::basic().unwrap();
        // base = 0.5 
        // Complete infinite sum
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.75);

        // [0,infinity]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 0), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.5);
        
        // [-infinity, 0]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity), &Float::with_val(32, 0)).unwrap();
        assert_eq!(s,0.5);

        // [1,infinity]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 1), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.25);

        // [-1,infinity]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, -1), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.625);

        // [-infinity,-1]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity),&Float::with_val(32, -1)).unwrap();
        assert_eq!(s,0.25);

        // [-infinity,1]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity),&Float::with_val(32, 1)).unwrap();
        assert_eq!(s,0.625);

        // [1,5]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 1),&Float::with_val(32, 5)).unwrap();
        assert_eq!(s,0.2421875);

        // [-5,5]
        let base = arithmeticconfig.get_float(0.5);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, -5),&Float::with_val(32, 5)).unwrap();
        assert_eq!(s,0.734375);


        // base = 0.1 
        // Complete infinite sum
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.55);

        // [0,infinity]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 0), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.5);
        
        // [-infinity, 0]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity), &Float::with_val(32, 0)).unwrap();
        assert_eq!(s,0.5);

        // [1,infinity]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 1), &Float::with_val(32, Special::Infinity)).unwrap();
        assert!(s.to_f64()-0.05 < 0.001); // result not exact, test case only

        // [-1,infinity]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, -1), &Float::with_val(32, Special::Infinity)).unwrap();
        assert_eq!(s,0.545);

        // [-infinity,-1]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity),&Float::with_val(32, -1)).unwrap();
        assert!(s.to_f64()-0.05 < 0.001); // result not exact, test case only

        // [-infinity,1]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, Special::NegInfinity),&Float::with_val(32, 1)).unwrap();
        assert_eq!(s,0.545);

        // [1,5]
        let base = arithmeticconfig.get_float(0.1);
        let s = get_sum(&base, &arithmeticconfig, &Float::with_val(32, 1),&Float::with_val(32, 5)).unwrap();
        assert!(s.to_f64()- 0.499995< 0.001);// result not exact, test case only


    }

    // Basic test of special infinity values
    #[test]
    fn test_infinity(){
        let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
        let a = arithmeticconfig.enter_exact_scope();
        assert!(a.is_ok());
        let f = arithmeticconfig.get_float(Special::Infinity);
        assert!(f.is_infinite());
        assert!(f.is_sign_positive());
        let g = arithmeticconfig.get_float(Special::NegInfinity);
        assert!(g.is_infinite());
        assert!(g.is_sign_negative());
        let b = arithmeticconfig.exit_exact_scope();
        assert!(b.is_ok());
        assert!(!g.is_integer()); // infinite values are not considered integers
    }

}
