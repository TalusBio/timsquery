
# Timsquery

UNDER CONSTRUCTION

... The idea is to have a library that allows querying TIMS data in a generic way.

The main design is to have modular aggregators, indices and queries.
Thus, depending on the access pattern, the data can be queried in
different ways and aggregated in different ways (if you need random
acces, use the index that works for that, if you need bulk
sequential, use that).

Ideally we will also have a sane CLI (think ... msaccess-like).
