#!/usr/bin/env python
from extractor import Extractor
from input import InputReader, CQPFormat
import utils

class SentenceLengthExtractor(Extractor):
    
    def getName(self):
        return "SentenceLengthExtractor"
    
    def getFeatureNames(self):
        return self.__featureNames
    
    def getAvailableOptions(self):
        short = ""
        long = []
        return short, long
        
    def __init__(self):
        self.__featureNames = utils.getStatNames("SENT_LENGTH")
        Extractor.__init__(self)
        
    def process(self,file):
        feats = {}
        Extractor.process(self,file)
        ir = InputReader(file)
        ir.read()
        cqpf = CQPFormat(ir.getText())
        lengths = [end-start for (start,end,arg) in cqpf.getAnnotations("s")]
        print self.__featureNames
        feats = utils.getStats("SENT_LENGTH", lengths)
        return ir.getID(),feats
        
if __name__ == "__main__":
    SentenceLengthExtractor()