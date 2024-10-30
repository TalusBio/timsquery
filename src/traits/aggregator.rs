/// A trait that defines how to aggregate items.
///
/// The `Item` type is the type of the item that is being aggregated.
/// The `Output` type is the type of the output of the aggregation.
///
/// The `add` method takes an item of type `Item` OR a type that
/// imlements `Into<Item>`.
///
/// The `finalize` method returns the output of the aggregation.
pub trait Aggregator: Send + Sync {
    type Item: Send + Sync + Clone;
    type Output: Send + Sync;

    fn add(&mut self, item: impl Into<Self::Item>);
    fn finalize(self) -> Self::Output;
}
