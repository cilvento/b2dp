use criterion::{criterion_group, Criterion, BenchmarkId};
use b2dp::{Eta, GeneratorOpenSSL,
    utilities::bounds::PartitionBound,
    utilities::bounds::PartitionBoundOptions,
    utilities::weights::WeightTable};




fn run_mechanism(x: &Vec<i64>, pb: &PartitionBound) -> u32 {
    let eta = Eta::new(1,1,1).unwrap();
    let _wt = WeightTable::from_bounds(eta, &pb, &x).unwrap();
    0
}


fn bench_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("LaplaceBoundWeightComputation");
    for n in [10,20,30,40,50].iter() {
        // Construct target partition of appropriate size
        let x: Vec<i64> = (0..*n).map(|x| x*x).rev().collect();
        let total: i64 = x.iter().sum();
        println!("{:?} sum: {:?}",x,total);
        let rng = GeneratorOpenSSL {};
        let eta = Eta::new(3,2,1).unwrap();
        let options: PartitionBoundOptions = Default::default();
        let pb = PartitionBound::from_noisy_estimates(total as usize, Some(*n as usize), &x, eta, rng, options).unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(total) , &total, 
            |b, n| b.iter(|| run_mechanism(&x, &pb)));
    }
    group.finish();
}

criterion_group!(benches, bench_sizes);
//criterion_main!(benches);