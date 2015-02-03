import os
import multiprocessing as mp

## see: http://matthewrocklin.com/blog/work/2013/12/05/Parallelism-and-Serialization/
## and https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma
## and https://stackoverflow.com/questions/19984152/what-can-multiprocessing-and-dill-do-together
##import pathos.multiprocessing as mp

import numpy as np
import pandas as pd
import sqlite3

print 'importing utils'

def stop():
    raise 'STOPPING ON PURPOSE!'

def system( cmd ):
    print cmd
    tmp = os.popen( cmd ).read()
    return tmp

def do_something_par( items, func, threads=None ): # None ensures that it makes as many threads as avail. in computer
    if threads == 1:
        out = map(func, items)
    else:
        pool = mp.Pool(processes=threads)              # start 4 worker processes
        ## Need to send, e.g. a tuple (1, counts_g) if fill_all_scores_par() took multiple args
        out = pool.map(func, items)
        pool.terminate()
    return out

def writeLines( lines, fname=None ):
    if fname is None:
        import tempfile
        fname = tempfile.mktemp()
    handle = open(fname, 'w')
    handle.write('\n'.join(lines)) ## assumes lines are an array split by '\n' - if not then do '\n'.join(lines) first
    handle.close()
    return fname

def readLines( fname ):
    fo = open( fname, 'r' )
    lines = fo.readlines()
    fo.close()
    return lines

def table( arr ):
    from collections import Counter
    c = Counter( arr )
    return c ## can be used like a dict.

def reverse_dict( d ): ## reverse keys <-> elements for a dict
    out = { v:k for k,v in d.items() }
    return out
 
def slice_sampler(px, N = 1, x = None):
    """
    Provides samples from a user-defined distribution.
    
    slice_sampler(px, N = 1, x = None)
    
    Inputs:
    px = A discrete probability distribution.
    N  = Number of samples to return, default is 1
    x  = Optional list/array of observation values to return, where prob(x) = px.
 
    Outputs:
    If x=None (default) or if len(x) != len(px), it will return an array of integers
    between 0 and len(px)-1. If x is supplied, it will return the
    samples from x according to the distribution px.    

    From:
    http://www.adamlaiacano.com/post/14987215771/python-function-for-sampling-from-an-arbitrary-discrete
    """
    from numpy.random import uniform
    import random
    values = np.zeros(N, dtype=np.int)
    samples = np.arange(len(px))
    px = np.array(px) / (1.*sum(px))
    u = uniform(0, max(px))
    for n in xrange(N):
        included = px>=u
        choice = random.sample(range(np.sum(included)), 1)[0]
        values[n] = samples[included][choice]
        u = uniform(0, px[included][choice])
    if x:
        if len(x) == len(px):
            x=np.array(x)
            values = x[values]
        else:
            print "px and x are different lengths. Returning index locations for px."
    if N == 1:
        return values[0]
    return values

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
        Taken from here: https://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

## TODO: look into bottleneck package for faster functions:
##   https://pypi.python.org/pypi/Bottleneck
def sdize_vector( vec, ignore_zeroes=True, use_median=True ): ## note this is inplace! If don't want, pass vec.copy() !!
    v = vec
    if ignore_zeroes:
        v = vec[ vec != 0 ]
    if use_median:
        from scipy.stats import nanmedian
        mn = nanmedian(v)
        sd = mad(v)
    else:
        mn = np.nanmean( v )
        sd = np.nanstd( v )
    vec -= mn
    vec /= (sd + 0.001) ## try to minimize copies?
    return vec

def minmax( vec ):
    return (np.nanmin(vec), np.nanmax(vec))

from matplotlib import pyplot as plt

def setup_text_plots(fontsize=8, usetex=True):
    """
    This function adjusts matplotlib settings so that all figures in the
    textbook have a uniform format and look.
    """
    import matplotlib
    matplotlib.rc('legend', fontsize=fontsize, handlelength=3)
    matplotlib.rc('axes', titlesize=fontsize)
    matplotlib.rc('axes', labelsize=fontsize)
    matplotlib.rc('xtick', labelsize=fontsize)
    matplotlib.rc('ytick', labelsize=fontsize)
    matplotlib.rc('text', usetex=usetex)
    matplotlib.rc('font', size=fontsize, family='sans-serif',
                  style='normal', variant='normal',
                  stretch='normal', weight='normal')

def read_all_tables( dbfile ):
    """read out all tables in a sqlite3 db file into a dict of pandas dataframes"""
    conn = sqlite3.connect( dbfile )
    tables = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name", conn)
    query = 'select * from %s'
    tables = { tname: pd.read_sql(query % tname, conn) for tname in tables.name.values }
    return tables
