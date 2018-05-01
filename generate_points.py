#!/usr/bin/env python
#
# generate_points.py
#
# Generate some integral points in a given bounding box and write them to
# stdout.
#

import sys
import time
import random

def usage(errstr=''):
    sys.stderr.write("""usage: %s [-s<seed>] [-w<width>=800 -h<height>=800] <N>

Generate N integral points in a given bounding box and write them to stdout.
Default width and height are 800. If either of width or height are given, both
must be given.  Default seed is epoch time in seconds.
""" % (sys.argv[0],))
    if errstr:
        sys.stderr.write('\nError: ' + errstr + '\n')
    sys.exit(1)

def generate_points(os, npoints, width=800, height=800):
    for pnum in xrange(npoints):
        px = random.randrange(width)
        py = random.randrange(height)
        os.write('{} {}\n'.format(px, py))

def option(optdict, sopt, lopt):
    if sopt in optdict:
        return optdict[sopt]
    if lopt in optdict:
        return optdict[lopt]
    return None

def main(*args):
    seed = None

    # Parse options and args
    from getopt import getopt, GetoptError
    sopts = 'hs:W:H:'
    lopts = ['help', 'seed=', 'width=', 'height=']
    o, args = getopt(args, sopts, lopts)
    o = dict(o)

    if option(o, '-h', '--help'):
        usage()

    if len(args) != 1:
        usage('Wrong number of arguments.')
    seed = int(option(o, '-s', '--seed') or time.time())
    width = int(option(o, '-W', '--width') or 800)
    height = int(option(o, '-H', '--height') or 800)
    npoints = int(args[0], 0)

    random.seed(seed)
    generate_points(sys.stdout, npoints, width, height)
    return 0

if __name__ == '__main__':
    errstr = ''
    ret = 0
    # Catch bad integer arguments
    try:
        ret = main(*sys.argv[1:])
    except ValueError as ve:
        errstr = str(ve)

    if errstr:
        sys.stderr.write('Error: ' + str(errstr) + '\n')
        if ret == 0:
            ret = 1
    sys.exit(ret)
