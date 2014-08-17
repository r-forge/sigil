"""featex_v2 -- Linguistic Feature Extraction Framework in Python Reloaded.

Preprocessor base class.

A preprocessor is similar to a feature extractor, but differs in the fact
that it is always called with ALL files of the corpus, and that it does 
not write its output in a predefined way.

"""

__version__ = "$Id: preprocessor.py 43 2007-06-30 15:54:55Z sm $"

import sys
from glob import glob
from getopt import *
from os.path import exists
import logging
import hashlib

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s: %(levelname)s %(message)s')

class Preprocessor:
    
    def __init__(self):
        #self.__name = "Extractor"
        self.logger = logging.getLogger(self.getName())
        self.__parseArgs()
        self.run()
        self.finish()
        
    def __parseArgs(self):
        short, long = self.getAvailableOptions()
        long.append('help')
        try:
            opts, args = getopt(sys.argv[1:],'h' + short,long)
        except GetoptError, e:
            return
        
        for (opt,arg) in opts:
            if opt=='-h' or opt=='--help':
                self.usage()
                return
        self.parseOpts(opts)
        files = []
        for a in args:
            files.extend(glob(a))
        self.files = [f for f in files if exists(f)]
        if len(self.files)==0:
            self.usage();
            print "Error: No files to process!"
            return
    
    def getAvailableOptions(self):
        return '',[]
    
    def parseOpts(self,opts):
        pass
    
    def finish(self):
        pass
    
    def getName(self):
        return "Preprocessor"
    
    def usage(self):
        print self.__doc__    
        
    def run(self):
        for file in self.files:
            self.process(file)
            
    def pr(self,text):
        self.logger.info(text)
            
    def process(self,file):
        self.pr("Processing file " + file)