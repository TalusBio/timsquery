use crate::traits::aggregator::{Aggregator, ProvidesContext};

pub trait QueriableData<QF, I, C>
where
    QF: Send + Sync + ProvidesContext<Context = C>,
    I: Send + Sync + Clone + Copy,
{
    fn query(&self, fragment_query: &QF) -> Vec<I>;
    fn add_query<A, O, AG>(&self, fragment_query: &QF, aggregator: &mut AG)
    where
        A: From<I> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O, Context = C>;
    fn add_query_multi_group<A, O, AG>(&self, fragment_queries: &[QF], aggregator: &mut [AG])
    where
        A: From<I> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O, Context = C>;
}
