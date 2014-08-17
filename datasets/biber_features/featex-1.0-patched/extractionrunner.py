#!/usr/bin/env python
"""featex_v2 -- Linguistic Feature Extraction Framework in Python Reloaded.

ExtractionRunner, run feature extraction with many extractors and many files.
"""

from getopt import *
from glob import glob
import sys
from os.path import exists, abspath
from os import popen
import logging
import threading

# extractors to be run (their filenames, without the .py)
EXTRACTORS = []
FILES = []
# number of worker threads
NUM_WORKERS = 2

class Worker(threading.Thread):
    
    def __init__(self, files, extractors):
        threading.Thread.__init__(self)
        self.files = files
        self.extractors = extractors
        self.logger = logging.getLogger("exRunner." + self.getName())
    
    def run(self):
        self.logger.debug("Starting Thread, processing " + str(self.files) + " with extractors " + str(self.extractors))
        self.process()
        self.logger.debug("Stopping Thread, processing " + str(self.files) + " with extractors " + str(self.extractors))
        
    def process(self):
        for extractor in self.extractors:
            prog = popen(abspath(extractor) + ' ' +' '.join(self.files))
            for line in prog:
                self.logger.debug("Extrator Output: " + line)

def parseArgs():
    try:
        opts, args = getopt(sys.argv[1:],'he:',['help','extractors=','logfile='])
    except GetoptError, e:
        return
    
    LOGFILE = False
    
    for (opt,arg) in opts:
        if opt=='-h' or opt=='--help':
            self.usage()
            return
        if opt=='-e' or opt=='--extractors':
            exts = arg.split(",")
            if len(exts) != 0:
                for ext in exts:
                    EXTRACTORS.append(ext.strip())
        if opt=='--logfile':
            LOGFILE = True
            logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s: %(levelname)s %(message)s',filename=arg,filemode='w')
    
    if not LOGFILE:
        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s: %(levelname)s %(message)s')
    
        
    files = []
    for a in args:
        files.extend(glob(a))
    files = [f for f in files if exists(f)]
    FILES.extend(files)
    if len(FILES)==0:
        usage();
        print "Error: No files to process!"
        return

# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/425397
def split_seq(seq, size):
    if size > len(seq):
        size = len(seq)
    newseq = []
    splitsize = 1.0/size*len(seq)
    for i in range(size):
        newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
    return newseq
   
def usage():
    print __doc__

def main():
    parseArgs()
    logger = logging.getLogger("exRunner")
    logger.info("Files: " + str(FILES))
    logger.info("Extractors: " + str(EXTRACTORS))
    splits = split_seq(FILES,NUM_WORKERS)
    print splits
    for i in range(len(splits)):
        w = Worker(splits[i],EXTRACTORS)
        w.start()
if __name__ == "__main__":
    main()