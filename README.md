# The Voronoi Game

This project aims to implement some algorithms in the field of
Competitive Facility Location (CFL). There are many papers on the
problem itself (see [References](#References)), but there appear to be
few implementations floating around.

# Table of contents

<!--ts-->
   * [The Voronoi Game](#the-voronoi-game)
   * [Table of contents](#table-of-contents)
      * [Background](#background)
      * [Approach](#approach)
      * [L1 P2 Solution](#l1-p2-solution)
      * [Time-based P2 solution](#time-based-p2-solution)
   * [Shapefiles](#shapefiles)
   * [Building](#building)
   * [Authors](#authors)
   * [References](#references)

<!-- Added by: foreese, at: 2019-01-09T19:06-0500 -->

<!--te-->

## Background

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
that the identical solution in the L2 norm runs in
*\&Omicron(n<sup>2</sup>)* time. This is not unexpected for a small
problem in 2D, but can add up quite quickly when the input consists of
millions of customers. Input of that size is not always necessary; but
the input size depends on the domain.

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

## Approach

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
we replaced the concepts of distance in the max-depth with time.

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
*isochrome*, or a *polyline* such that every point on the boundary of
the isochrome takes a fixed time *t* to reach from the center. Each
round of the game, the max-depth algorithm will request (1) the travel
time from each customer to its nearest facility, and (2) an isochrome
matching that travel time centered at each customer. To produce these we
linearly interpolate between the known fixed rings according to the
Euclidean distance between the facility point and its nearest two rings,
or between the smallest ring and the center point. Past the end of all
known rings we extrapolate similarly.

It is possible to increase the accuracy of the algorithm by writing a
tool to query travel time offline, for example using an (“OpenStreetMap”
[2019](#ref-OpenStreetMap)) database or ArcGIS. Another improvement
would be to write a script to prefetch actual data from (“TravelTime
Platform” [2019](#ref-TravelTimePlatform)) – for now, the only data
tested with is randomly generated by a helper script in
`[generate.py](/scripts/generate.py)`

# Shapefiles

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

# Building

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

# Authors

**Fritz Reese**  
*M.S. Computer Science*, George Mason University, 2018.

**[Jyh-Ming Lien](https://cs.gmu.edu/~jmlien/doku.php)**  
*Ph.D. Computer Science*, Texas A\&M University, 2006.  
Department of Computer Science  
George Mason University

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
