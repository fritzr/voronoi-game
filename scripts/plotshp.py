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

def plot_reader(reader):
    if reader.shapeType in (shapefile.POINT, shapefile.POINTZ):
        plot_points(reader, marker='o', markersize=4)
    elif reader.shapeType in (shapefile.POLYGON, shapefile.POLYGONZ):
        plot_polygons(reader)
    else:
        raise RuntimeError('unsupported shape type '+str(reader.shapeType))

def main(*args):
    p = argparse.ArgumentParser(description='Plot a shapefile.')
    p.add_argument('shapefiles', metavar='PATH...', type=str, nargs='+',
            help='Basenames to shapefiles.'
                 ' *.{shp,shx,dbf} should exist for each PATH.')

    paths = p.parse_args(args).shapefiles
    excs = list()
    pyplot.figure()

    for path in paths:
        try:
            reader = shapefile.Reader(path)
            plot_reader(reader)
        except Exception as e:
            excs.append(e)
        finally:
            if int(shapefile.__version__.split('.')[0]) > 1:
                reader.close()

    ret = 0
    if excs:
        ret = 1
        sys.stderr.write('Ignored the following errors:\n')
        for exc in excs:
            sys.stderr.write('  %s: %s\n' % (type(exc).__name__, str(exc)))

    pyplot.show()

    return ret

if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
