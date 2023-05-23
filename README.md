# Extended Morton Coder

A rust library that computes extended morton codes (see
[the abstract](https://dl.acm.org/doi/10.1145/3105762.3105782)). Instead of just
using the XYZ coordinates of the centroid, it additionally uses the length of
the diagonal to split off primitives with different sizes. It also has adaptive
axis order, to use more bits in a skewed dimension.

Reason to use size bits - if you can split off large primitives from the smaller
ones, you can keep them together in a BVH, away from the other primitives. This
means the large AABB of the large primitives will not "infect" the smaller ones,
allowing you to ignore a large number of small primitives with fewer AABB tests.

Reason to use adapative axis order - if your scene is skewed, say it's 4x longer
in the X and Y dimensions than in the Z dimension, you can allocate more bits
for the X and Y part of the code, enabling better tracing performance.

# Notes

Still a WIP. See issues for various bugs. This crate is in a more-or-less
finished state (although PRs are still welcome to fix those issues!)
