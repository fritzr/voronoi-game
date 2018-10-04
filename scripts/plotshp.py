#!/usr/bin/env python

import sys
import argparse
import shapefile
from matplotlib import pyplot

def plot_points(reader, *pltargs, **pltkws):
    for shape in reader.shapes():
        x, y = shape.points[0]
        pyplot.plot(x, y, *pltargs, **pltkws)

def plot_polygons(reader, *pltargs, **pltkws):
    for polygon in reader.shapes():
        x = [pt[0] for pt in polygon.points]
        y = [pt[1] for pt in polygon.points]
        pyplot.plot(x, y, *pltargs, **pltkws)

def main(*args):
    p = argparse.ArgumentParser(description='Plot a shapefile.')
    p.add_argument('shapefile', metavar='PATH', type=str,
            help='Basename to shapefile. Expects PATH.{shp,shx,dbf} to exist.')

    path = p.parse_args(args).shapefile
    reader = shapefile.Reader(path)
    try:
        pyplot.figure()
        if reader.shapeType in (shapefile.POINT, shapefile.POINTZ):
            plot_points(reader, marker='o', markersize=4)
        elif reader.shapeType in (shapefile.POLYGON, shapefile.POLYGONZ):
            plot_polygons(reader)
        else:
            raise RuntimeError('unsupported shape type '+str(reader.shapeType))
        pyplot.show()
    except:
        raise
    finally:
        if int(shapefile.__version__.split('.')[0]) > 1:
            reader.close()

    return 0

if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
