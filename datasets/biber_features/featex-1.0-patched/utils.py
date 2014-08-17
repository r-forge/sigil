from getopt import *
import sys
import math

def parseOpts(shortopts,longopts):
        """Parse options with getopt. """
        try:
            opts, args = getopt(sys.argv[1:],shortopts,longopts)
        except GetoptError, e:
            return
        return opts, args
    
def median(list):
    """Calculate the median value of the given list of numbers"""
    if len(list) % 2 == 0:
        s = sorted(list)
        return 0.5 * s[len(list)/2 - 1] +  0.5 * s[len(list)/2]
    else: 
        return sorted(list)[len(list)/2]
    
def interQuartileRange(list):
    s = sorted(list)
    return s[3 * len(list)/4] - s[len(list)/4]  
    
def mean(list):
    return sum(list) / float(len(list))

def variance(list):
    s = sum(list)
    sumsq = sum([i**2 for i in list])
    return ((sumsq - s**2/float(len(list)))/float(len(list)))

def powsum(list,power):
    return sum([i**power for i in list])

def moment(list,mean,power):
    return sum([(i - mean)**power for i in list])

def stdev(list):
    return math.sqrt(variance(list))

def getStatNames(varname):
    return [varname + "_" + st for st in STATS.keys()]

def getStats(varname, sample):
    out = {}
    for s in STATS:
        out[varname + "_" + s] = STATS[s](sample)
    return out

STATS = {"mean":mean, "median":median, "min":min, "max":max, "var":variance, "stdev": stdev, "iquart":interQuartileRange}
