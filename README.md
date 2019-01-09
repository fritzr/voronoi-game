# Competitive Facility Location with the Voronoi Game

This project aims to implement some algorithms in the field of
*Competitive Facility Location*. There are many papers on the problem
itself (see [References](#References)), but there appear to be few
implementations floating around.

# Table of contents

<!--ts-->

  - [Competitive Facility Location with the Voronoi
    Game](#competitive-facility-location-with-the-voronoi-game)
  - [Table of contents](#table-of-contents)
  - [Input](#input)
  - [Building](#building)
  - [Authors](#authors)
  - [References](#references)

<!-- Added by: fritz, at: 2018-10-15T18:35-0400 -->

<!--te-->

# Methods

This project is an attempt at implementing a practical and realistic
*Competitive Facility Location* algorithm. It is based on several works,
primarily those of (Imai and Asano 1983) and (Banik et al. 2017).

# Input

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
<https://doi.org/https://doi.org/10.1016/j.comgeo.2017.02.003>.

</div>

<div id="ref-IMAI1983310">

Imai, Hiroshi, and Takao Asano. 1983. “Finding the Connected Components
and a Maximum Clique of an Intersection Graph of Rectangles in the
Plane.” *Journal of Algorithms* 4 (4): 310–23.
<https://doi.org/https://doi.org/10.1016/0196-6774(83)90012-3>.

</div>

</div>
