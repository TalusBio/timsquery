use crate::Aggregator;

pub trait QueriableData<QF, I>
where
    QF: Send + Sync,
    I: Send + Sync + Clone + Copy,
{
    fn query(&self, fragment_query: &QF) -> Vec<I>;
    fn add_query<A, O, AG>(&self, fragment_query: &QF, aggregator: &mut AG)
    where
        A: From<I> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O>;
    fn add_query_multi_group<A, O, AG>(&self, fragment_queries: &[QF], aggregator: &mut [AG])
    where
        A: From<I> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O>;
}

// I like this idea but I need a way to set/propagate the tolerance
//
// impl<T: ToleranceAdapter<QF, ElutionGroup<K>>, QF, A, K> QueriableData<QF, A> for T {
//     fn query(&self, fragment_query: &QF) -> Vec<A> {
//         let mut out = Vec::new();
//         let qf = self.query_from_elution_group(fragment_query);
//         self.add_query(&qf, &mut |x| out.push(x));
//         out
//     }
//
//     fn add_query<O, AG: Aggregator<Item = A, Output = O>>(
//         &self,
//         fragment_query: &QF,
//         aggregator: &mut AG,
//     ) {
//         let qf = self.query_from_elution_group(fragment_query);
//         self.add_query(&qf, aggregator);
//     }
//
//     fn add_query_multi_group<O, AG: Aggregator<Item = A, Output = O>>(
//         &self,
//         fragment_queries: &[QF],
//         aggregator: &mut [AG],
//     ) {
//         let qfs = fragment_queries
//             .iter()
//             .map(|x| self.query_from_elution_group(x))
//             .collect();
//         self.add_query_multi_group(&qfs, aggregator);
//     }
// }
