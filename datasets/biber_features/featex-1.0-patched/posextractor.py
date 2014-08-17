#!/usr/bin/env python
from extractor import Extractor
from input import InputReader, CQPFormat
import utils

class POSExtractor(Extractor):
    """Extractor for POS uni-, bi-, and trigrams.

    Usage: posextractor.py -f [POSFILE] [FILES]

    [POSFILE] has to contain one uni-, bi-, or trigram per line, with
    single tags separated by underscores, e.g. "AT0", "AT0_NN0", "AT0_AJ0_NN1".
    
    The [FILES] have to contain one token per line, with columns separated
    by a tab-character. The first column has to contain the word, the 
    second column the POS-tag (from the BNC C5 tagset), and the third
    column the lemmatized form.

    The output is written to a file named POSExtractor_[UNIQUEID].tab,
    where [UNIQUEID] is determined from the file names given. This allows
    several extraction proccesses to be started on subsets of the data
    without interference.
    """    
    def getName(self):
        return "POSExtractor"
    
    def getFeatureNames(self):
        return self.__featureNames
    
    def getAvailableOptions(self):
        short = "f:"
        long = ["posfile="]
        return short, long
    
    def parseOpts(self,opts):
        self.posfile = "posextract.txt"
        for (opt,arg) in opts:
            if opt == "-f" or opt == "--posfile":
                self.posfile = arg
        self.readNames()
        self.__featureNames = self.posnames
                
    
    def readNames(self):
        print self.posfile
        f = open(self.posfile)
        self.posnames = [l.strip() for l in f]
        f.close()
        
    def __init__(self):
        Extractor.__init__(self)
        
    def process(self,file):
        feats = {}
        Extractor.process(self,file)
        ir = InputReader(file)
        ir.read()
        cqpf = CQPFormat(ir.getText())
        pos = cqpf.getColumn(1)
        # initialize counts
        
        for name in self.posnames:
            feats[name] = 0
        
        for i in range(2,len(pos)): # ignore first two pos ...
            uni =  (pos[i])[0:3]
            bi = (pos[i-1])[0:3] + "_" + uni
            tri = (pos[i-2])[0:3] + "_" + bi
            if uni in feats:
                feats[uni] += 1
            if bi in feats:
                feats[bi] += 1
            if tri in feats:
                feats[tri] += 1
            
        for x in self.posnames:
            feats[x] /= float(len(pos)-2)
        
        return ir.getID(),feats
        
if __name__ == "__main__":
    POSExtractor()
