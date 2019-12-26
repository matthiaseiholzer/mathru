
#[cfg(test)]
mod normal
{
    use mathru::statistics::distrib::{Continuous, Normal};

    #[test]
    fn pdf0()
    {
        let mean : f64 = 0.0;
        let variance: f64 = 1.0;
        let distrib : Normal = Normal::new(mean, variance);
        let x : f64 = 0.0;
        let prob : f64 = distrib.pdf(x);

        assert_eq!(0.3989422804014327, prob);
    }

// Does not work all the time, because the used function random is not mocked.
//    #[test]
//    fn random()
//    {
//        let mean_1 : f64 = 0.0;
//        let variance_1: f64 = 1.0;
//        let distrib_1 : NormalDistrib = NormalDistrib::new(&mean_1, &variance_1);
//        let mut data: Vec<f64> = Vec::new();
//
//        for _i in 0..10000
//        {
//            data.push(distrib_1.random());
//        }
//
//        let distrib_2: NormalDistrib = NormalDistrib::from_data(&data);
//        println!("{}", distrib_2.variance());
//        assert!(distrib_2.mean() < 0.01);
//        assert!(distrib_2.mean() > -0.01);
//        assert!(distrib_2.variance() < 1.02 * variance_1);
//        assert!(distrib_2.variance() > 0.98 * variance_1);
//    }

    #[test]
    fn cdf0()
    {
        let mean: f64 = 0.0;
        let variance: f64 = 1.0;
        let distrib: Normal = Normal::new(mean, variance);

        assert_eq!(0.5, distrib.cdf(0.0))
    }

    #[test]
    fn quantile0()
    {
        let mean: f64 = 0.0;
        let variance: f64 = 1.0;
        let distrib: Normal = Normal::new(mean, variance);

        assert_eq!(1.2815515655446006, distrib.quantile(0.9));
    }

    #[test]
    fn quantile1()
    {
        let mean: f64 = 1.0;
        let variance: f64 = 0.5;
        let distrib: Normal = Normal::new(mean, variance);

        assert_eq!(1.0, distrib.quantile(0.5));
    }

//    #[test]
//    fn from_data()
//    {
//        let mean: f64 = 5.0;
//        let variance: f64 = 10.0;
//        let num_samples: usize = 100;
//        let data: Vector<f64> = Normal::new(mean, variance).random_vector(num_samples);
//
//        let distrib: Normal = Normal::from_data(&data);
//
//        assert!((mean - distrib.mean()).abs() < 0.5);
//        assert!((variance - distrib.variance()) < 1.0);
//    }
}