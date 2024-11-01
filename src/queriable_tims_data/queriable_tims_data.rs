use rayon::prelude::*;
use serde::Serialize;
use std::hash::Hash;
use std::time::Instant;
use tracing::info;

use crate::traits::aggregator::ProvidesContext;
use crate::{Aggregator, ElutionGroup, QueriableData, Tolerance, ToleranceAdapter};

// TODO: URGENTLY make documentation fot eh functions using this leftover struct docs
// as a reference.
//
// A struct that can be queried for TIMS data.
//
// The main idea behind this struct is to provide a generic way to query TIMS data.
// And by that I mean that with the same interface but different imeplementations
// on the indexing one can optimize for different access patterns.
//
// In the same way different aggregators can be used to aggregate the results.
// In some cases just adding the intensities might be enough, but in other cases
// one might want to do more complex operations. (fitting a gaussian to a chromatogram for example)
//
// Main Generic parameters:
// - `QD`: The type of indexed data that will be queried.
//         This can also optimize the access pattern to the data. When the data is queried
//         with the `many` version of the query functions, the indexed data can optimize
//         the access pattern to the data.
//
//         also implements tolerance adapter that will be used to convert elution groups
//         into queries.
// - `TL`: The type of tolerance that will be used to define the search space.
//
// Additional Generic parameters:
// - `QP`: The type of precursor query that will be used to query the indexed data.
// - `QF`: The type of fragment query that will be used to query the indexed data.
// - `OE`: The type of output element that will be returned by the aggregators.
// - `AG`: The type of aggregator that will be used to aggregate the output elements.
// - `AE`: The type of element that will be aggregated.
//
//
// So ... in other words ...
// 1. `TA` converts `ElutionGroup` with `TL` into `QP` and `QF` queries.
// 2. `QD` is queried with `QP` and `QF` queries.
// 3. The results are aggregated with `AG` aggregators (which are generated by the factory, when passed a numeric ID).

pub fn query_multi_group<'a, QD, TL, QF, AE1, AE2, OE, AG, FH, CTX1, CTX2, FF>(
    queriable_data: &'a QD,
    tolerance: &'a TL,
    elution_groups: &[ElutionGroup<FH>],
    aggregator_factory: &FF,
) -> Vec<OE>
where
    AG: Aggregator<Item = AE2, Output = OE, Context = CTX2> + Send + Sync,
    QD: QueriableData<QF, AE1, CTX1> + ToleranceAdapter<QF, ElutionGroup<FH>>,
    TL: Tolerance,
    OE: Send + Sync,
    FH: Clone + Eq + Serialize + Hash + Send + Sync,
    QF: Send + Sync + ProvidesContext<Context = CTX1>,
    // AE: Send + Sync + Clone + Copy,
    AE1: Into<AE2> + Send + Sync + Clone + Copy,
    AE2: Send + Sync + Clone + Copy + From<AE1>,
    CTX1: Into<CTX2> + Send + Sync + Clone + Copy,
    FF: Fn(u64) -> AG + Send + Sync,
{
    let start = Instant::now();
    let mut fragment_queries = Vec::with_capacity(elution_groups.len());
    let mut aggregators = Vec::with_capacity(elution_groups.len());

    for (i, elution_group) in elution_groups.iter().enumerate() {
        let qp = queriable_data.query_from_elution_group(tolerance, elution_group);

        fragment_queries.push(qp);
        aggregators.push(aggregator_factory(i as u64));
    }

    queriable_data.add_query_multi_group(&fragment_queries, &mut aggregators);
    let duration = start.elapsed();
    info!("Querying took {:#?}", duration);

    let microsecond_duration = duration.as_micros();
    let microseconds_per_query = microsecond_duration / elution_groups.len() as u128;
    let queries_per_second = 1_000_000.0 / microseconds_per_query as f64;

    info!("That is {:#?} queries per second", queries_per_second);
    info!(
        "That is {:#?} microseconds per query",
        microseconds_per_query
    );

    let start = Instant::now();
    let out = aggregators.into_par_iter().map(|x| x.finalize()).collect();
    let elapsed = start.elapsed();
    info!("Aggregation took {:#?}", elapsed);
    out
}
