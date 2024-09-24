pub trait Aggregator<I>: Send + Sync {
    type Output: Send + Sync;

    fn add(&mut self, item: &I);
    // fn fold(&mut self, item: Self);
    fn finalize(self) -> Self::Output;
}
