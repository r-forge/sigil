"""featex_v2 -- Linguistic Feature Extraction Framework in Python Reloaded.

Feature Extractor base class.

"""

__version__ = "$Id: extractor.py 80 2007-07-27 19:13:53Z sm $"

import sys
from glob import glob
from getopt import *
from os.path import exists
from time import asctime
import logging
import hashlib

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s: %(levelname)s %(message)s')

class Extractor:
    
    def __init__(self):
        #self.__name = "Extractor"
        self.logger = logging.getLogger(self.getName())
        self.feats = {}
        self.__parseArgs()
        self.featfile = self.getFeatFile()
        self.writeHeader()
        self.run()
        
    def __parseArgs(self):
        short, long = self.getAvailableOptions()
        long.append('help')
        try:
            opts, args = getopt(sys.argv[1:],'h'+short,long)
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
        
    def getFeatFile(self):
        m = hashlib.md5()
        m.update(''.join(self.files))
        fn = self.getName() + "_" + m.hexdigest() + ".tab"
        f = open(fn,"w")
        return f
    
    def writeHeader(self):
        self.featfile.write('#ID\t')
        self.featfile.write('\t'.join(self.getFeatureNames()))
        self.featfile.write('\n')
    
    def writeFeats(self,id,featDict):
        out = []
        out.append(id)
        for feat in self.getFeatureNames():
            out.append(repr(featDict[feat]))
        self.featfile.write('\t'.join(out))
        self.featfile.write('\n')
        self.featfile.flush()
    
    def getName(self):
        return "Extractor"
    
    def getFeatureNames(self):
        return []
    
    def usage(self):
        print self.__doc__    
        
    def run(self):
        for file in self.files:
            id, feats = self.process(file)
            self.feats[id] = feats
            self.writeFeats(id,feats)
    
    def pr(self,text):
        self.logger.info(text)
            
    def process(self,file):
        self.pr("Processing file " + file)