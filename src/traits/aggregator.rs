pub trait Aggregator: Send + Sync {
    type Item: Send + Sync;
    type Output: Send + Sync;

    fn add(&mut self, item: &Self::Item);
    fn finalize(self) -> Self::Output;
}
