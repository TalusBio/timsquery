pub trait Aggregator<I>: Send + Sync {
    type Output;

    fn add(&mut self, item: &I);
    // fn fold(&mut self, item: Self);
    fn finalize(self) -> Self::Output;
}
