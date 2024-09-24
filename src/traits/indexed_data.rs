use crate::Aggregator;

pub trait IndexedData<QF, A> {
    fn query(&self, fragment_query: &QF) -> Vec<A>;
    fn add_query<O, AG: Aggregator<A, Output = O>>(&self, fragment_query: &QF, aggregator: &mut AG);
    fn add_query_multi_group<O, AG: Aggregator<A, Output = O>>(
        &self,
        fragment_queries: &[QF],
        aggregator: &mut [AG],
    );
}
