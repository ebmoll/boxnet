boxnet
======

a 2D sorting spatial subdivision algorithm

Boxnet is written in C and has a very simple to use API.
It can be used as a drop-in replacement for other broad-phase
algorithms like BSP, Quadtree, spatial hash, sort-and sweep or
BB-Tree.

Boxnet does not use a hierarchical subdivision, nor does it use
fixed-size cells. Instead it cuts the 2D space into non-overlapping,
variable-sized rectangles, and stores the nearest neighbors for each
rectangle for traversal.


see "BUILDING" for building information
