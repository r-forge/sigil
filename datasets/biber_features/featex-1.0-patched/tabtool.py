#!/usr/bin/env python
"""tabtool.py -- tool for handling tab separated files with headers in the
first line and IDs in the first column.

Usage: tabtool.py [operation] [file1] [file2] ... [fileN] [result]

Supported operations are: 

    cat   --  concatenate the rows of the given files, keeping just one header
    merge --  merge the columns of the given files on the first column
    sort  --  sort file by first column
    auto  --  do cat and merge based on file name patterns
"""

from glob import glob
import sys
from os.path import basename,exists
from shutil import copyfile

def main():
    mode,fileArgs,outfile = parseArgs()
    files = []
    for a in fileArgs:
        files.extend(glob(a))
    files = [f for f in files if exists(f)]
    if len(files) < 2:
        print __doc__
        print "ERROR: Need at least 2 files to work on!"
        exit()
    if mode == "cat":
        cat(files,outfile)
    if mode == "merge":
        merge(files,outfile)
    if mode == "auto":
        auto(files,outfile)
    
def cat(files,outfile):
    print "Concatenating " + prettyprint(files) + " --> " + outfile
    out = open(outfile,"w")
    haveWrittenHeader = False
    header = ""
    for file in files:
        f = open(file)
        oldheader = header
        header = f.readline()
        if not haveWrittenHeader:
            out.write(header)
            haveWrittenHeader = True
        else:
            if header != oldheader:
                raise ValueError(file + ": Headers do not match!")
        try:
            for line in f:
                out.write(line)
        finally:
            f.close()
    out.close()

def merge(files,outfile):
    print "Merging " + prettyprint(files) + " --> " + outfile
    # Assuming entries occur in the same order!
    # Also assuming IDs are exactly 3 characters long!
    out = open(outfile,"w")
    fs = [open(f,"r") for f in files]
    line = fs[0].readline()[0:-1]
    while line != "":
        out.write(line)
        line = fs[0].readline()[0:-1]
        for i in range(1,len(fs)):
            out.write('\t')
            out.write(fs[i].readline()[4:-1])
            out.write('\n')
    for f in fs:
        f.close()
    out.close()

def prettyprint(files):
    out = ['\n']
    for file in files:
        out.append('\t')
        out.append(file)
        out.append('\n')
    return ''.join(out)    
    
def copy(infile,outfile):
    print "Copying " + infile + " --> " + outfile
    copyfile(infile,outfile)

def auto(files,outfile):
    print "Automagic processing:  " + prettyprint(files) + " --> " + outfile
    tocat = {}
    for file in files:
        name = basename(file).split('_')[0]
        if name not in tocat:
            tocat[name] = []
        tocat[name].append(file)
    for k in tocat:
        if len(tocat[k]) > 1:
            cat(tocat[k],k+".tab")
        else:
            copy(tocat[k][0],k+".tab")
    tomerge = []
    for k in tocat:
        tomerge.append(k+".tab")
    if len(tomerge) > 1:
        merge(tomerge,outfile)
    else:
        copy(tomerge[0],outfile)

def parseArgs():
    if len(sys.argv) < 4:
        print __doc__
        print "ERROR: Need at least 3 arguments!"
        exit()
    else:
        mode = sys.argv[1]
        if mode not in ['cat', 'merge', 'auto']:
            print __doc__
            print "ERROR: Mode " + mode + " not recognized!"
            exit()
        return (mode,sys.argv[2:-1],sys.argv[-1])

if __name__ == "__main__":
    main()

