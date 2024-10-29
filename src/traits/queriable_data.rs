use super::tolerance::ToleranceAdapter;
use crate::models::elution_group::ElutionGroup;
use crate::Aggregator;

pub trait QueriableData<QF, A>
where
    QF: Send + Sync,
    A: Send + Sync + Clone + Copy,
{
    fn query(&self, fragment_query: &QF) -> Vec<A>;
    fn add_query<O, AG: Aggregator<Item = A, Output = O>>(
        &self,
        fragment_query: &QF,
        aggregator: &mut AG,
    );
    fn add_query_multi_group<O, AG: Aggregator<Item = A, Output = O>>(
        &self,
        fragment_queries: &[QF],
        aggregator: &mut [AG],
    );
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
