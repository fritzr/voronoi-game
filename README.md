# The Voronoi Game

This problem aims to implement a practical algorithm in the field of
Competitive Facility Location (CFL) in *R<sup>2</sup>*, also known as
*The Voronoi Game*.

# Table of Contents

<!--ts-->
   * [The Voronoi Game](#the-voronoi-game)
   * [Table of Contents](#table-of-contents)
   * [Background](#background)
      * [Competitive Facility Location](#competitive-facility-location)
      * [Example](#example)
      * [CFL in <em>R<sup>2</sup></em>](#cfl-in-r2)
      * [Line Segment Intersection](#line-segment-intersection)
   * [Approach](#approach)
      * [L1 P2 Solution](#l1-p2-solution)
      * [Time-based P2 solution](#time-based-p2-solution)
         * [Plane Sweep](#plane-sweep)
            * [Plane Sweep Issues](#plane-sweep-issues)
            * [Plane Sweep Runtime](#plane-sweep-runtime)
         * [Max Depth Calculation](#max-depth-calculation)
            * [Max Depth Issues](#max-depth-issues)
   * [Future Work](#future-work)
   * [Building](#building)
      * [CMake](#cmake)
      * [Dependencies](#dependencies)
      * [Introspection and Debug mode](#introspection-and-debug-mode)
   * [Running](#running)
      * [Shapefiles](#shapefiles)
      * [CLI](#cli)
   * [Authors](#authors)
   * [References](#references)

<!-- Added by: foreese, at: 2019-01-13T16:05-05:00 -->

<!--te-->

# Background

## Competitive Facility Location

The CFL problem, also called the Voronoi Game, is to optimally
distribute resources (facilities) assigned to different entities
(players) given a customer space and behavioral model. In this game,
each player takes turns with rounds of placing facilities in the
customer space in order to service the most number of customers
possible. The problem becomes a geometric one when the criteria for
selecting a facility is related to distance. Often, we assume that a
customer will attend the closest facility, since it is most convenient.

The problem has a rich history in various domains, particularly in the
transportation and GIS world. An excellent review of the literature in
its many domains is presented by (Eiselt, Laporte, and Thisse
[1993](#ref-EISELT199344)).

## Example

A good example domain is pizza delivery: consider two pizza companies,
Pizza John and Papa Hut, both of whom are competing to build stores in a
geographic region with the goal of maximzing profit. Assuming both
companies offer pizza of equal quality, with equal service, and at the
same times, it is almost certain that the customers in the region will
order pizza from whichever store is closest, to obtain their pizza
soonest. At the start of the game, Pizza John may already have 4 stores
in the region when Papa Hut seeks to enter the market. The game starts
when the “second player”, Papa Hut, places a new store. Obviously Papa
Hut will try to place one or more stores such that the number of
customers closer to Papa Hut stores than Pizza John stores is maximized.

## CFL in *R<sup>2</sup>*

If we consider the customer space in *R<sup>2</sup>* as in (Banik et al.
[2017](#ref-BANIK201753)), where the customer distance metric is
straight-line Euclidean distance, the answer is clearly related to the
concept of Voronoi Diagrams. Consider a Voronoi Diagram such that the
center of the cells (sites) are set to the location of the facilities
(stores). Each Voronoi cell can be thought of as “owned by” the player
who owns the facility at its center. The player’s “score” is defined as
the number of customers lying inside the player’s Voronoi cells
constructed with the player’s facilities at the sites.

(Banik et al. [2017](#ref-BANIK201753)) pursues this analogy and
describes the obvious solution in *R<sup>2</sup>*. Consider drawing a
circle centered on a customer point such that the radius of the circle
is the distance from that customer to its nearest facility. By
definition, any other point that lies on the interior of the circle will
be closer to the customer than the original facility. Therefore, by
constructing this map with circles centered about each customer point,
the ideal region for a new facility is the region which is formed by the
interior of the most number of intersecting circles. We call this the
“max-depth intersecting region”, and it is the solution in
*R<sup>2</sup>* for player 2.

Before proposing the *R<sup>2</sup>* solution using circles and
straight-line Euclidean distance, (Banik et al.
[2017](#ref-BANIK201753)) describes a simpler approximation method which
uses Manhattan distance. In the L1 norm, the fixed-distance boundary
shape is a tilted square rather than a circle. Banik shows one can
rotate the coordinate space then easily compute the max-depth region on
the resulting axis-aligned squares, resulting in an identically-aligned
rectangle. Banik refers to (Imai and Asano [1983](#ref-IMAI1983310)) for
this method. The general idea is to record intersections of rectangles
in an adjacency list using a plane sweep algorithm. Then, compute the
max clique of the resulting “connected components” graph. In general, a
graph may have up to *3<sup>n/3</sup>* maximal cliques as described in
(De, Nandy, and Roy [2014](#ref-DE201458)); but Imai & Asano show that
due to the special characteristics of a graph of intersecting
rectangles, the solution is always obtainable in *Ο(n* log *n)* time.

This runtime is excellent, but the solution is approximate. Banik shows
that the identical solution in the L2 norm runs in *Ο(n<sup>2</sup>)*
time. This is not unexpected for a small problem in 2D, but can add up
quite quickly when the input consists of millions of customers. Input of
that size is not always necessary; but the input size depends on the
domain.

A related and more complex problem is the player 1 solution. That is,
where should Pizza John, who is already in the market, place its first 4
stores to minimize the maximal gain from a competitor who enters the
market, assuming the competitor will follow the player strategy outlined
above? When there are no competitor facilities, the ideal location for
the first facility is the planar center (and its extensions) described
by (Matousek [1991](#ref-MATOUSEK2000221)). However, it becomes more
complicated when both players already have facilities in the region.
Banik offers a solution for this problem in *R<sup>2</sup>*, using a
similar circle method, which runs in *Ο(n<sup>8</sup>)* time.

Others propose methods outside *R<sup>2</sup>* ((Teramoto, Demaine, and
Uehara [2011](#ref-JGAA-235)), (Banik et al. [2016](#ref-BANIK201641)),
(Drezner [2014](#ref-DREZNER2014)), (Eiselt and Laporte
[1989](#ref-EISELT1989231))).

## Line Segment Intersection

As we will see in [Approach](#approach), we will use a modified version
of the classical plane sweep algorithm for computing line segment
intersections. Here we review the classical plane sweep approach for
reference.

The input is a set of line segments. We maintain the event queue *Q*,
the sweep status *T*, and the sweep line *L*.

*Q* is a priority queue containing events, which are either endpoint
events (endpoints of the input segments) or intersection events (points
where segments intersect). The queue is initialized with all endpoints
of the input segments, and intersection points are discovered as *L*
sweeps the plane. Event points are visited in order of descending *y*
coordinate. At each event we draw the sweep line *L* at the same *y*
value as the event point. *T* is a binary search tree that stores
references to all segments which intersect *L* sorted by the *x*
coordinate at which each segment intersects *L*. Note that the event
points are sorted lexicographically, so all points at the same *y*
coordinate are visited in order of increasing *x* coordinate.

The fundamental observation the drives the algorithm is that for
segments to intersect, their *y* intervals and *x* intervals have to
overlap. By restricting the status tree *T* to segments that intersect
*L*, we know their *y* intervals overlap. Then by sorting *T* by *x*
coordinate of all segments at each *L*, we know immediately that only
adjacent segments in *T* can possibly intersect. The trick is that, at
each intersection point formed by two segments *A* and *B*, the
positions of *A* and *B* in *T* are swapped. At this moment we can check
again whether *A* and *B* intersect the next adjacent edges in *T* and
queue any further intersection points. To maintain the invariant that
*T* holds every segment that intersects *L*, we insert a segment to *T*
when we visit its top endpoint and remove the segment from *T* when we
visit its bottom endpoint.

# Approach

This project aims to provide a CFL implementation with realistic results
and constraints. Though the *Ο(n<sup>8</sup>)* runtime of the player 1
solution is theoretically “reasonable” (it is polynomial, afterall) the
bound is likely too steep to run on a large input space (millions of
customers) in a reasonable amount of time. Furthermore, some
implementation details were left unclear from (Banik et al.
[2017](#ref-BANIK201753)). For these reasons and the sake of time we
decided to implement the simpler player 2 (P2) strategy of max-depth
intersection (with a few modifications).

## L1 P2 Solution

First, we implemented the player 2 solution in the L1 norm (Manhattan
distance) as described in (Imai and Asano [1983](#ref-IMAI1983310)). We
used a range-tree approach to find each customers’ nearest facility
based on (Cabello et al. [2010](#ref-CABELLO201099)), then followed Imai
& Asano’s paper using a line sweep. In the end, we avoided the special
tree structure *T* from Imai & Asano to track the depth of the
intermediate regions, as it would require manual manipulation of nodes
in a custom 2-3 tree. A reliable 2-3 tree implementation with direct
access to the internal nodes could not be found, and we preferred to
spend more time on a more realistic solution, so we avoided creating our
own 2-3 tree structure. Therefore we opted for a more input-sensitive
technique which can obtain *Ο(n* log *n)* in the average case but
*Ο(n<sup>2</sup>)* in the worst case, where we manually increment the
depths of inner rectangles when an encompassing rectangle is inserted.
The code for this is found in `src/maxrect.{cc,h}` from the `MaxRect`
class.

The results of this method are mixed. For small enough input sets it can
provide a reasonably accurate solution, but for larger data sets the
solution can be significantly worse than ideal.

After implementing the L1 solution with mixed results (the L1 solution
can be significantly off from the L2 solution) we decided to take a step
back and look at a more realistic method.

## Time-based P2 solution

First, we reconsidered the idea of straight-line distance as a customer
choice model. We realized that a distance metric is actually just an
approximation of *time*, which most directly represents the customers’
selection criteria (think again of the pizza example). To realize this,
we replaced the concepts of distance in the max-depth algorithm with
time.

There are two main requirements which change significantly when
replacing distance with time: first we need to identify the nearest
facility *in time* for every customer; second we need to identify a
shape which defines a fixed-*time* boundary isomorphic to the circles of
L2 and squares of L1.

Both data are actually readily available for real-world geographic
locations through services such as (“OpenStreetMap”
[2019](#ref-OpenStreetMap)) and (“TravelTime Platform”
[2019](#ref-TravelTimePlatform)). The latter supports exactly the
queries described above, but only through an online API. Since we want
to avoid web queries during active runs of a CPU-intensive geometric
algorithm, we decided to take as input a quantized cache of this data
for the query points. Each customer is input to the program as a
combination point and list of a fixed number of rings. Each ring is an
[isochrone](https://wiki.openstreetmap.org/wiki/Isochrone): in our usage
a polygon such that every point on the boundary of the isochrone takes a
fixed time *t* to reach from the center point. Each round of the game,
the max-depth algorithm will request (1) the travel time from each
customer to its nearest facility, and (2) an isochrone matching that
travel time centered at each customer. To produce these we linearly
interpolate between the known fixed rings according to the Euclidean
distance between the facility point and its nearest two rings, or
between the smallest ring and the center point. Past the end of all
known rings we extrapolate similarly.

See the section on [Future Work](#future-work) for a discussion on this
input method.

Following the (Imai and Asano [1983](#ref-IMAI1983310)) algorithm, the
algorithm can be broken into two main parts: the [plane
sweep](#plane-sweep) and [max depth
calculation](#max-depth-calculation).

The idea is that the maximal-depth region formed by the most number of
intersecting polygons will be a solution isomorphic to the maximal
intersecting number of squares or circles. Instead of intersecting
polygons, we decided to triangulate the polygons and intersect triangles
instead. The idea is that the arrangement of triangles (and the
resulting regions) is easier to calculate and represent than the
arrangement of polygons.

Unlike Banik’s solution with circles, we believe we can compute the
max-depth player 2 solution using this triangulated isochrone method in
an average time proportional to *n* log *n* (for *n* customers).

### Plane Sweep

First we need to detect the intersections of the triangles, as we would
with squares or circles. But the algorithm to detect intersections must
be slightly different when using isochrones instead of circles or
rectangles.

First, we require any continuous isochrones to be discretized into
polygons. Then, to simplify the intersection algorithm, we decompose the
polygons into triangles. Triangles intersect iff their edges intersect,
so we consider the edge set of the trianges as line segments and compute
the intersections of all segments.

In our first pass implementation to find the segment intersection points
we follow the classical [line segment
intersection](#line-segment-intersection) plane sweep algorithm, where
the event queue *Q* is a `std::priority_queue` of points, the sweep
status *T* is a `std::set` of segments, and the sweep line *L* is stored
as the last *y* coordinate visited. We use a `std::set` because its
search, removal, and insertion operations have logarithmic complexity
bounds. (In practice the `std::set` is usually a red-black tree).

To properly maintain the segments in *T*, a custom comparator is used
which intersects each segment with *L* and sorts by the *x*-value of the
intersection. At intersection points, both segments then intersect *L*
at the same *x* value, so we sort by orientation such that the left-most
oriented segment comes first. To see why, imagine that *L* moves an
infinitesimally small distance below the intersection point. Clearly the
left-oriented segment will intersect *L* at a lower *x* value.

At endpoint events, we update *L* and then insert/remove the segment
(for top/bottom endpoints). At intersection events, to swap segments we
remove the segments from *T*, update *L*, then re-insert the segments
(with a location hint for efficiency) so they are ordered correctly past
the intersection point.

We accomplish this using `std::set` by providing the comparator with an
indirect reference to *L*, so the sorting criteria changes when *L* is
updated. Strictly speaking this breaks the invariant that the comparator
object and elements of a `std::set` remain `const`. However, the only
time the order of segments in *T* can ever change is at an intersection
point, due to the design of the algorithm. At this point the only
segments whose order in *T* could be affected is the segments that make
up the intersection. Since we remove them before updating *L*, *T*
remains consistent.

#### Plane Sweep Issues

It turns out that this implementation is not robust for intersection
points which are formed by more than two segments. Unfortunately, for
our input, this is guaranteed to happen, because our line segments are
from triangles which form triangulated polygons. Every triangle will
share edges with two other triangles; that is, we will have many
segments which appear twice (once for each triangle) but are really the
same edge/segment. We can’t completely collapse the edge into a single
segment because we need to acknowledge the intersection of any lines
through both triangles separately for the [max depth
calculation](#max-depth-calculation).

Consider if we have a set of segments *S* which all intersect at point
*P*. To handle this, we need more than a single “swap”. The most
intuitive way to think about this is that each pair from
*℘<sub>2</sub>(S)* (all pairs possible from *S*) forms a separate
intersection event consisting of a single swap. This follows the
operation of the classical algorithm, BUT the pairwise swaps must occur
in the correct order; that is, in order of orientation from left to
right. Then each segment in *S* must be checked for intersections with
the adjacent segments in *T* that are *not* in *S* (line segments cannot
intersect each other twice).

So far it appears this intuitive approach *could* work. However,
consider the case where two identical segments *B* and *C* are
intersected by some non-colliner segment *A*, and that the order of
these segments is in the status tree *T* is *(A, B, C)*. In theory, if
*B* is ordered before *C* for whatever reason, first the intersection
event *A ∩ B* is queued. When the intersection event is visited, *A* and
*B* are swapped so *T* contains *(B, A, C)*. At this point, *A* and *B*
are checked for intersections with their neighbors, and *A ∩ C* is found
and queued. The intersection is the next event handled (since it’s at
the same point) at which time *A* will swap with *C* and *T* will
contain *(B, C, A)*.

In order to implement the algorithm so the above example works, since we
simply remove and re-insert segments to perform swaps, at some point the
comparator would have to decide to sort *A* between *B* and *C*. But *B*
and *C* are identical segments, which intersect *L* at the same *x*
value and have the same orientation\! Furthermore, *A*’s orientation is
to the right of both *B* and *C*, so the comparator could never decide
that *T* should contain *(B, A, C)*. It would immediately sort *A* after
*C* after the first removal and re-insertion when handling *A ∩ B*. In
this case, *C* is never checked for intersections with adjacent segments
in the status. The problem is exacerbated if there are additional
segments which intersect at the same point.

To solve this, it is clear that all segments that form an intersection
point must be “handled” in a single visitation of the intersection
event. That is, when a set of segments *S* all intersect at point *P*,
we should handle *P* once and in so doing check every segment in *S* for
intersections with the left- and right-adjacent segments in *T* which
are *not* in *S*.

To that end, our implementation maintains an additional structure *I*
which is an associative container (`std::unordered_map`) mapping
intersection points *P* to the set *S<sub>P</sub>* of segments
(`std::unordered_set`) which intersect at *P*. The order of segments in
the set is not important, since the segments themselves are still sorted
properly in *T*. When it is discovered that two segments *X* and *Y*
intersect at a point *P*, the point is hashed, used to find
*S<sub>P</sub>* from *I*, and then *X* and *Y* are inserted into
*S<sub>P</sub>*. Since it is a set, no segment will appear twice. When
an intersection event *E* is handled, the corresponding set
*S<sub>E</sub>* is looked up from *I*. All segments in *S<sub>E</sub>*
are located in and removed from *T* before *L* is updated. Once *L* is
updated, all segments are inserted again (in no particular order) and
checked for intersections with the segments in *T* adjacent but outside
*S<sub>E</sub>*. This way, duplicates are collapsed and no segments are
skipped in the re-sorting step.

#### Plane Sweep Runtime

It is worth discussing the runtime of the sweep algorithm after our
modifications. Let *n* be the number of segments in the input. Clearly
|*T*| is bounded by *n*. The cost for an endpoint event is one insertion
(or removal) in *T* each requiring time *Ο(* log *n)* followed by one or
two checks for intersection with the neighboring segments in *T*.
Accessing a neighboring segment and checking two segments for
intersection take amortized constant time. Therefore the total cost for
all endpoint events is *Ο(n* log *n)*.

The cost for an intersection event *E* is one removal and one insertion
in *T* and one or two segment intersection checks for each segment in
*S<sub>E</sub>*. The insertion and removals again take time *Ο(* log
*n)* and the adjacent traversals and intersection checks are amortized
constant time. It follows that the time for each intersection event is
*Ο(k<sub>E</sub>* log *t)*, where *k<sub>E</sub> = *|*S<sub>E</sub>*|.
Let the total time for all intersection events be *T<sub>I</sub>*. Let
*K =* Σ*<sub>E</sub> k<sub>E</sub>*. Then *T<sub>I</sub> ≤ K* log *n*.

Let *m* be the number of intersection events. *m* is upper-bounded by
Σ<sub>*i=1,n*</sub> *i ≤ M* where *M = n<sup>2</sup>/2*. To see this,
consider adding segments to the input one-by-one so that each insertion
creates the maximal number of intersections; the *i+1*-th segment can
contribute *i* intersections by intersecting all the *i* existing
segments. Any inserted segment that forms an intersection with other
segments at an existing intersection point *P* instead adds one to *K*
and contributes only *i - k<sub>P</sub>* to *m*. By summing the terms we
see that *m ≤ M - K* and *m + K ≤ M*. Thus, the amortized traversal over
*k<sub>E</sub>* segments for each intersection event *E* does not cost
us anything more than the classical traversal of only pairwise
intersection events. Therefore the runtime of our sweep algorithm is
still *Ο((n + m)* log *n)* which is input-sensitive and bounded by
*Ο(n<sup>2</sup>* log *n)* as with the classical plane sweep algorithm
for line segment intersections.

### Max Depth Calculation

Once the intersections between all triangulated customer isochrones are
computed, we can go ahead and compute the maximal depth region. Our
original approach for this was to follow (Imai and Asano
[1983](#ref-IMAI1983310)) again. In their paper, visiting all the
intersections of the input shapes (rectangles) computed the connected
components graph. This graph is represented as an adjacency list (or
matrix) wherein each vertex corresponds to an input shape, and an edge
between two vertices indicates that the two corresponding shapes
intersect (overlap). In the resulting graph, every maximally connected
subgraph represents a single “connected component”. Imai & Asano prove
that maximal cliques in this graph correspond 1:1 with the maximal depth
intersecting regions formed by the input rectangles.

To copy this procedure for triangles, we set up an adjacency matrix
where the indexes correspond to the indexes of the input triangles.
During the [plane sweep](#plane-sweep) algorithm, whenever we handle an
intersection between segments we simply look up the owning triangles and
mark an edge between their corresponding vertices in the adjacency
matrix.

The idea is that the maximal clique for such an intersection graph of
triangles, like the corresponding graph of rectangles, represents
triangles which intersect to form a maximal-depth region.

#### Max Depth Issues

Unfortunately, after some thought, it becomes clear that the solution
Imai & Asano present for rectangles does not translate cleanly to
triangles. It is true that any maximal-depth region with depth *d* will
form a clique of degree *d*. However, it is easy to imagine *d+1*
triangles which all intersect each other, but do not intersect in any
single closed region with depth *d*. In fact you can construct any
number of triangles which all intersect each other such that no more
than two triangles intersect in the same region. So even by searching
among maximal cliques of depth *d+1* or greater one will never find the
true maximal-depth region of depth *d*.

# Future Work

The first step moving forward is to rewrite the [max depth
calculation](#max-depth-calculation) algorithm. Unfortunately, as
discussed above, the max clique technique will not work for triangles.
The correct algorithm would be to compute the entire arrangement of the
plane: that is, extend the plane sweep algorithm to build the faces of
all sub-polygons formed by the intersections of the triangles or
polygons, keeping track of the “depth” of each region. Faces must be
split and merged as necessary, reconciling the depth at each step. One
must handle intersections by forming a new polygonal face with depth
combined from the two intersecting edges’ shapes. With this method it is
not immediately clear whether the breakdown into triangles really
simplifies the algorithm (from a runtime complexity standpoint). At
first thought it certainly seems simpler – fewer cases to handle of
intersecting edges. It is also worth noting that the max-depth region
will always be convex, and it may be that not *every* face needs be
stored throughout the entire sweep. For the results of this project to
be correct, this algorithm must be thoroughly designed and implemented
correctly.

An additional improvement beyond this would be to obtain and process
real input data. As mentioned above, (“OpenStreetMap”
[2019](#ref-OpenStreetMap)) and (“TravelTime Platform”
[2019](#ref-TravelTimePlatform)) are excellent resources for such data.
Instead of obtaining arbitrary discrete isochrones and interpolating
between them, the ideal scenario is to write a library to query the
exact travel time isochrones offline.

On improvement to note is that this project takes a shortcut when
“generating” an interpolated isochrone. At each round, every user
should have an isochrone representing the fixed travel time to its
nearest facility. This is isomoprhic to the fixed-distance squares or
circles from the L1 and L2 solutions. Though our program correctly
interpolates between isochrones to approximate the travel time, it does
not generate a new isochrone at that time. Such an isochrone
interpolation algorithm was outside the scope of this project, but is
feasible and there are discussions on the web for how to implement them.
ArcGIS itself has the capability for this, and there is a Python project
called [isochrones](https://github.com/timothydmorton/isochrones) which
may provide useful results.

# Building

## CMake

This project builds with CMake to help locate dependencies. This means
usually the following is
sufficient:

``` bash
$ mkdir build && cd build && cmake [flags] .. && make -j4 && make -j4 install
```

If you want to control the build, you may need to modify some cmake
build parameters. I find the CMake documentation incredibly dense, so
here are a few common cmake flags for reference:

``` 
  -DCMAKE_INSTALL_PREFIX= # like --prefix for autoconf
  -DCMAKE_CXX_FLAGS=      # build flags passed to CXX
```

## Dependencies

Version 2.8.12 is the minimum tested CMake version, but \>= 3.0 should
be fine too.

The source code uses C++11 features. This means `g++ >= 4.8.5` and a
functional C++ build environment is required. If you are using `g++`
versions 5.0 through 5.4, you will also need static C++ libraries (on
linux, it is typically sufficient to install `libstdc++-static`.

Additionally, development versions of the following libraries are
required to build. The versions shown are the earliest known accepted
versions:

| Package  | Known Supported Versions |
| -------- | ------------------------ |
| OpenCV   | \>= 2.4.5                |
| boost    | \>= 1.59.0               |
| GNU GMP  | \>= 6.0.0                |
| GNU MPFR | \>= 3.1.1                |
| Shapelib | \>= 1.3.0                |

To use the shapefile helpers in the `scripts/` subdirectory you will
need:

| Package | Known Supported Versions |
| ------- | ------------------------ |
| Python  | 2.7 - 3.x                |
| pyshp   | \>= 1.2.1                |

Additionally, for `plotshp.py`, you will need `tkinter` and
`matplotlib`.

This project also uses a library called [FIST: Fast Industrial-Strength
Triangulation of
Polygons](http://www.cosy.sbg.ac.at/~held/projects/triang/triang.html)
(Copyright Martin Held) for triangulating the polygons. The source code
is not in the public domain, therefore only the headers and a static
library binary are made available in this project (see
[include/FIST/](include/FIST) and [lib/FIST/](lib/FIST)). Unfortunately
the FIST library I have access to is built for 32-bit x86, which may
limit the distributability of this project. However, if you need to
rebuild FIST, you may obtain the source code for academic or commercial
use by [contacting Martin
Held](https://www.cosy.sbg.ac.at/~held/email.html).

## Introspection and Debug mode

There are a few project-specific flags you can pass to cmake:

``` 
  -DDEBUG=[0/1]
  -DMAXTRI_DEBUG=[0/1]
```

If you pass `-DDEBUG=1`, extra debugging output will be enabled.

If you pass `-DMAXTRI_DEBUG=1`, a LOT of extra debugging output will be
enabled for the max-depth triangle algorithm. Included in this is the
dumping of the status tree during the plane sweep to Octave/MATLAB
script files in the current directory for each step of the status tree.

# Running

## Shapefiles

This program reads and writes in the [ESRI
Shapefile](https://en.wikipedia.org/wiki/Shapefile) format. This format
is chosen because it is widely used in the GIS community, and the
specification is open. ESRI has published a
[whitepaper](https://www.esri.com/library/whitepapers/pdfs/shapefile.pdf)
detailing the file specification.

There are several ways to read and write these files. This project
provides some very basic Python scripts in the `scripts/` subdirectory
which can be used to generate and plot shapefile test data.

There are low-level APIs available for many languages; see for example
[shapelib](http://shapelib.maptools.org/) and
[GDAL](https://www.gdal.org/) for C/C++, or
[pyshp](https://github.com/GeospatialPython/pyshp) for Python.

The easiest and most powerful way to view and create shapefiles
(especially large files) is to use the proprietary
[ArcGIS](https://arcgis.com) software. If you don’t have access to
ArcGIS, you can use [GNU Octave](https://www.gnu.org/software/octave/)
with the
[octave-mapping](https://octave.sourceforge.io/mapping/index.html)
plugin package. (Note this package requires `octave-geometry` and
`octave-io`.)

A “shapefile” is actually a collection of several files with the same
base name and different extensions. The shape file itself ends in
`.shp`, but alongside it is usually at least `.dbf` and `.shx` files.
Whenever a “path to a shapefile” is requested, the path requested should
actually the basename common to the `.shp`, `.dbf` and `.shx` files.

The format of the input shapefiles is documented in the help text for
the [`generate.py`](scripts/generate.py) script. This script is used to
generate simple test input for the program:

    $ ./scripts/generate.py --help
    usage: generate.py [OPTIONS...] [-p|-u [POINT_FILE]] [-P [POLY_FILE]] N
    
    Output polygons centered around N random points to one or more Shapefiles.
    -p/-u writes points, -P writes polygons. If neither are given, write both
    to './points.*' and './polygons.*'.
    -p writes pure points (no fields), -u writes 'user' points (see below).
    
    Polygons are generated by drawing points on the boundary of an inner ring,
    then protruding a number of 'spike' points onto the boundary of an outer ring.
    The result is a star-like polygon.
    
    Output is written in the ESRI Shapefile format, see:
        http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
    
    With -u the point Shapefile describes users, containing the following records:
        name           type        value
        ================================
        pointIndex     N (double)  index of point
        weight         N (double)  normalized weight (0-1)
    
    The user points need not exactly represent individual persons; they should
    probably be compressed to represent some "population unit" in a raster-like
    field, where the weight could represent population density near the point.
    
    The polygon Shapefile contains the following records:
        name           type        value
        ================================
        pointIndex     N (int)     index of associated point in point file
        FTT            N (double)  fixed travel time from point to ring in minutes
    
    positional arguments:
      npoints               Number of points to generate.
    
    optional arguments:
      -h, --help            show this help message and exit
    
      Common options
    
      -s SEED, --seed SEED  Seed the RNG with the given integer. Default current
                            time.
      -x MIN_X, --min-x MIN_X
                            Min X coordinate of any point. Default 0.0.
      -X MAX_X, --max-x MAX_X
                            Max X coordinate of any point. Default 1920.0.
      -y MIN_Y, --min-y MIN_Y
                            Min Y coordinate of any point. Default 0.0.
      -Y MAX_Y, --max-y MAX_Y
                            Max X coordinate of any point. Default 1080.0.
      -o FILE, --output FILE
                            Write using basename FILE. Default is './shape'. Note
                            that the Shapefile format specifies the creation of
                            FILE.shp, FILE.shx, and FILE.dbf.
      -D, --debug           Debug with PDB.
      -p [FILE], --point-file [FILE]
                            Generate points to FILE.{{shp,shx,dbf}}.
      -P [FILE], --poly-file [FILE]
                            Generate polygons to FILE.{{shp,shx,dbf}}.
      -g, --graph           Graph points and polygons when writing them.
    
    Polygon parameters:
      -d DIRECTIONS, --directions DIRECTIONS
                            Number of directional spikes. Default is 5.
      -n VERTICES, --vertices VERTICES
                            Average number of vertices per polygon. Default is
                            twice the number of directions (from -d).
      -r INNER_RADIUS, --inner-radius INNER_RADIUS
                            Average radius of the inner ring. Default 40.0.
      -R OUTER_RADIUS, --outer-radius OUTER_RADIUS
                            Average radius of the outer ring. Default 80.0.
      -l LAYERS, --layers LAYERS
                            Number of 'layers' for each polygon. N layers of
                            concentric polygons are generated around each point.
                            Default 1.
      -L LAYER_FACTOR, --layer_factor LAYER_FACTOR
                            Additive scale factor between layers. Default 2.0.
      -T FTT_FACTOR, --ftt_factor FTT_FACTOR
                            Fixed Travel Time (FTT) between polygon layers, in
                            minutes. Default 1.0 minute(s).
      -t FTT, --ftt FTT     Initial Fixed Travel Time (FTT) for the first isoline,
                            in minutes. Default 1.0 minute(s).
      -v VARIANCE, --variance VARIANCE
                            Variance of polygon parameters. Parameter names should
                            be the same as the option names above in this section
                            (short or long) with dashes replaced by underscores.
                            The value given is the distance from the average to
                            the upper and lower boundary in the random number. (A
                            variance of 2 and an average of 5 is the random range
                            [3,7]). Default string is:
                            'a:35.0,d:2,n:3,r:5.0,R:10.0,l:0,L:0.0,T:0.0,t:0.0'
    
    In this example:
    
        {prog}  --variance=d:2,r:10,R:10,n:5  -d5 -r40 -R80 -n10  15
    
    we generate 15 polygons, each with their parameters randomly assigned in
    the range shown:
    
        inner radius: [30, 50]
        outer radius: [70, 90]
        polygon vertices: [5, 15]
        directional spikes: [3, 7]

## CLI

The program itself has a simple command-line interface. In its current
state, it is expected that extra code need be written to extend this
project. The command-line interface as it stands can be observed by
running the program with
    `--help`:

    usage: voronoi-game [OPTIONS] [-p<N>] <user_points> <user_rings> <p1_sites> <p2_sites> <output>
    
    Read users and their fixed-travel-time (FTT) rings, along with some known facilities/sites for players 1 and 2. Then successively find the optimal location to place sites for the second player. With -p, play N rounds, considering each player alternatively as 'player 2'. The output is a Point Shapefile containing the solution points.
    
    OPTIONS are:
      -D            Debug mode (extra output)
      -p N          Play N rounds of the game (default: 1).
      -s SEED       Seed the RNG with the given number.
                Default is to use current Epoch time.

As described above in the [Shapefiles](#shapefiles) section, you can see
the program expects the customer points and rings, and the player one
and two site files. The output file will be written to a simple Point
shapefile with no database entries (though a vacuously empty DBF file
will be written).

# Authors

  **Fritz Reese**  
  * *M.S. Computer Science*, George Mason University, 2018.

  **[Jyh-Ming Lien](https://cs.gmu.edu/~jmlien/)**  
  * [Motion and Shape Computing (MASC) group](http://masc.cs.gmu.edu/)
  * Associate Prof., [Department of Computer Science](https://cs.gmu.edu/), George Mason University

# References

<div id="refs" class="references">

<div id="ref-BANIK201753">

Banik, Aritra, Bhaswar B. Bhattacharya, Sandip Das, and Satyaki
Mukherjee. 2017. “The Discrete Voronoi Game in R2.” *Computational
Geometry* 63 (Supplement C): 53–62.
<https://doi.org/10.1016/j.comgeo.2017.02.003>.

</div>

<div id="ref-BANIK201641">

Banik, Aritra, Jean-Lou De Carufel, Anil Maheshwari, and Michiel Smid.
2016. “Discrete Voronoi Games and ϵ-Nets, in Two and Three Dimensions.”
*Computational Geometry* 55 (Supplement C): 41–58.
<https://doi.org/https://doi.org/10.1016/j.comgeo.2016.02.002>.

</div>

<div id="ref-CABELLO201099">

Cabello, S., J.M. Díaz-Báñez, S. Langerman, C. Seara, and I. Ventura.
2010. “Facility Location Problems in the Plane Based on Reverse Nearest
Neighbor Queries.” *European Journal of Operational Research* 202 (1):
99–106. <https://doi.org/https://doi.org/10.1016/j.ejor.2009.04.021>.

</div>

<div id="ref-DE201458">

De, Minati, Subhas C. Nandy, and Sasanka Roy. 2014. “In-Place Algorithms
for Computing a Largest Clique in Geometric Intersection Graphs.”
*Discrete Applied Mathematics* 178 (Supplement C): 58–70.
<https://doi.org/https://doi.org/10.1016/j.dam.2014.06.025>.

</div>

<div id="ref-DREZNER2014">

Drezner, Tammy. 2014. “Competitive Facility Location in the Plane” 7
(December).

</div>

<div id="ref-EISELT1989231">

Eiselt, H.A., and G. Laporte. 1989. “Competitive Spatial Models.”
*European Journal of Operational Research* 39 (3): 231–42.
<https://doi.org/https://doi.org/10.1016/0377-2217(89)90161-6>.

</div>

<div id="ref-EISELT199344">

Eiselt, H. A., Gilbert Laporte, and Jacques-François Thisse. 1993.
“Competitive Location Models: A Framework and Bibliography.”
*Transportation Science* 27 (1): 44–54.
<https://doi.org/10.1287/trsc.27.1.44>.

</div>

<div id="ref-IMAI1983310">

Imai, Hiroshi, and Takao Asano. 1983. “Finding the Connected Components
and a Maximum Clique of an Intersection Graph of Rectangles in the
Plane.” *Journal of Algorithms* 4 (4): 310–23.
<https://doi.org/10.1016/0196-6774(83)90012-3>.

</div>

<div id="ref-MATOUSEK2000221">

Matousek, Jiri. 1991. “Computing the Center of Planar Point Sets.”
*Papers from the DIMACS Special Year, Comput. Geom.*, 221–30.

</div>

<div id="ref-OpenStreetMap">

“OpenStreetMap.” 2019. <https://wiki.openstreetmap.org>.

</div>

<div id="ref-JGAA-235">

Teramoto, Sachio, Erik D. Demaine, and Ryuhei Uehara. 2011. “The Voronoi
Game on Graphs and Its Complexity.” *Journal of Graph Algorithms and
Applications* 15 (4): 485–501. <https://doi.org/10.7155/jgaa.00235>.

</div>

<div id="ref-TravelTimePlatform">

“TravelTime Platform.” 2019. <https://app.traveltimeplatform.com>.

</div>

</div>
