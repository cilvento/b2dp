use criterion::{criterion_group, Criterion, BenchmarkId};
use rug::{ Float, float::Special};
use b2dp::utilities::{exactarithmetic::ArithmeticConfig,discretesampling::get_sum};


fn run_sums(offset: u32,sum_type: u32) {

    let mut arithmeticconfig = ArithmeticConfig::basic().unwrap();
    arithmeticconfig.increase_precision(16).unwrap();
    
    match sum_type {
        1 => { // Complete infinite sum
            let base = arithmeticconfig.get_float(0.5);
            let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, Special::NegInfinity), &Float::with_val(arithmeticconfig.precision, Special::Infinity)).unwrap();
        },
        // 2 => { // [-infinity, 0]
        //     let base = arithmeticconfig.get_float(0.5);
        //     let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, Special::NegInfinity), &Float::with_val(arithmeticconfig.precision, 0.0)).unwrap();
        // }
        3 => { // [-infinity, offset]
            let base = arithmeticconfig.get_float(0.5);
            let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, Special::NegInfinity), &Float::with_val(arithmeticconfig.precision, offset)).unwrap();
        },
        4 => { // [-infinity, -offset]
            let base = arithmeticconfig.get_float(0.5);
            let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, Special::NegInfinity), &Float::with_val(arithmeticconfig.precision, -(offset as i64))).unwrap();
        },
        // 5 => { // [0,infinity]
        //     let base = arithmeticconfig.get_float(0.5);
        //     let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, 0), &Float::with_val(arithmeticconfig.precision, Special::Infinity)).unwrap();
        // },
        6 => { // [-offset,infinity]
            let base = arithmeticconfig.get_float(0.5);
            let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, -(offset as i64)), &Float::with_val(arithmeticconfig.precision, Special::Infinity)).unwrap();
        },
        7 => { // [offset,infinity]
            let base = arithmeticconfig.get_float(0.5);
            let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, offset), &Float::with_val(arithmeticconfig.precision, Special::Infinity)).unwrap();
        },
        8 => { // [-5,5]
            let base = arithmeticconfig.get_float(0.5);
            let _s = get_sum(&base, &arithmeticconfig, &Float::with_val(arithmeticconfig.precision, -(offset as i64)), &Float::with_val(arithmeticconfig.precision, offset + 1)).unwrap();
        },
        _ => {}
    };
}



fn bench_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("GetSum Comparisons");
    
    for i in 1..10 {
        group.bench_with_input(BenchmarkId::new("[-infinity, infinity]", i) , &i, 
            |b, i| b.iter(|| run_sums(*i,1)));
        group.bench_with_input(BenchmarkId::new("[-infinity, offset]", i) , &i, 
            |b, i| b.iter(|| run_sums(*i,3)));
        group.bench_with_input(BenchmarkId::new("[-infinity, -offset]", i) , &i, |b, i| b.iter(|| run_sums(*i,4)));
        group.bench_with_input(BenchmarkId::new("[-offset, infinity]", i) , &i, 
            |b, i| b.iter(|| run_sums(*i,6)));
        group.bench_with_input(BenchmarkId::new("[offset, infinity]", i) , &i, 
            |b, i| b.iter(|| run_sums(*i,7)));
        group.bench_with_input(BenchmarkId::new("[-offset, offset + 1]", i) , &i, 
            |b, i| b.iter(|| run_sums(*i,8)));
    }
    group.finish();
}



criterion_group!(benches, bench_sizes);
//criterion_main!(benches);