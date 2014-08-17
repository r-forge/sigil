#!/usr/bin/env python

from preprocessor import Preprocessor
from input import *
import utils
import psyco
psyco.full()
psyco.log()
psyco.profile(0.1)

class LexiconPP(Preprocessor):
    
    def getName(self):
        return "LexiconPreProcessor"
    
    def getAvailableOptions(self):
        return 'c:',['columns=']
    
    def parseOpts(self,opts):
        for (opt,arg) in opts:
            if opt == '-c' or opt == '--column':
                self.column = int(arg)
    
    def __init__(self):
        self.lexicon = {}
        self.count = 0
        self.column = 2
        Preprocessor.__init__(self)
    
    def process(self,file):
        Preprocessor.process(self,file)
        ir = InputReader(file)
        ir.read()
        cqpf = CQPFormat(ir.getText())
        for word in cqpf.getColumn(self.column):
            if word not in self.lexicon:
                self.lexicon[word] = 0
            self.lexicon[word] += 1
            self.count += 1
            
    def finish(self):
        f = open("lexicon.txt","w")
        for word in sorted(self.lexicon.keys()):
            f.write(word + "\t" + str(self.lexicon[word]) + "\t" + str(self.lexicon[word] / float(self.count)) + "\n")
            
if __name__ == "__main__":
    LexiconPP()
