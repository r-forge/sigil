#!/usr/bin/env python
# -*- coding: utf8 -*-
"""featex_v2 -- Linguistic Feature Extraction Framework in Python Reloaded

select.py -- Select (and symlink) documents based on certain meta-information 

Usage: select [-d datadir] [-m metafile] [-o outdir] expression
"""

import sys
import os, os.path
from getopt import getopt, GetoptError

## DEFAULT OPTION VALUES
DATADIR = "./"
METAFILE = "overview.tbl"
OUTDIR = "select_out/"
EXPRESSION ='True'

def usage():
    print __doc__

def lines2indexedDict(lines,idxcolumn,sep='\t'):
    out = {}
    heads = lines.next().split(sep)
    heads = [h.strip() for h in heads]
    numCols = len(heads)
    idxid = heads.index(idxcolumn)
    for line in lines:
        data = line.split(sep)
        if len(data) != numCols:
            raise ValueError("Number of columns not equal for head/data!")
        for i in range(0,numCols):
            if not data[idxid] in out:
                out[data[idxid]] = {}
            out[data[idxid]][heads[i]] = data[i]
    return out


def main():
    try:
        opts, args = getopt(sys.argv[1:],'hd:m:o:',['help','datadir=','metafile=','outdir='])
    except GetoptError, e:
        usage()
        print "Error: " + str(e)
        exit()
    if len(args) == 0:
        print "No expression given!"
        usage()
        exit()
    EXPRESSION = args[0]
    for (opt,arg) in opts:
        if opt=='-h' or opt=='--help':
            usage()
            return
        if opt=='-d' or opt=='--datadir':
            DATADIR = arg
        if opt=='-m' or opt=='--metafile':
            METAFILE = arg
        if opt=='-o' or opt=='--outdir':
            OUTDIR = arg
    print OUTDIR
    mf = open(METAFILE,"r")
    meta = lines2indexedDict(mf,'ID')
    if not os.path.isdir(DATADIR):
        print DATADIR + " is not an existing directory!"
        return
    DATADIR = os.path.abspath(DATADIR) + '/'
    if not os.path.isdir(OUTDIR):
        os.mkdir(OUTDIR)
    files = os.listdir(DATADIR)
    files = [f for f in files if os.path.isfile(DATADIR + f)]
    for f in files:
        fn = os.path.basename(f)
        dotpos = fn.find('.')
        id = fn[:dotpos]
        if not id in meta:
            print "Meta information for file " + f + " not found!"
            continue
        d = meta[id]
        if eval(EXPRESSION)==True:
            os.symlink(DATADIR + f,OUTDIR + f)
            print DATADIR + f + " -> " + OUTDIR + f
        
if __name__ == '__main__':
    main()
    
