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
