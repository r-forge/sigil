#!/usr/bin/env python

from preprocessor import Preprocessor
from input import *
import utils

class POSPP(Preprocessor):
    
    def getName(self):
        return "POSPreProcessor"
    
    def getAvailableOptions(self):
        return 'c:',['columns=']
    
    def parseOpts(self,opts):
        for (opt,arg) in opts:
            if opt == '-c' or opt == '--column':
                self.column = int(arg)
    
    def __init__(self):
        self.unilexicon = {}
        self.bilexicon = {}
        self.trilexicon = {}
        
        self.count = 0
        self.column = 1
        Preprocessor.__init__(self)
    
    def dis(self,pos):
        return pos[0:3]
    
    def process(self,file):
        Preprocessor.process(self,file)
        ir = InputReader(file)
        ir.read()
        cqpf = CQPFormat(ir.getText())
        pos = cqpf.getColumn(self.column)
        for i in range(2,len(pos)): # ignore first two pos ...
            uni =  (pos[i])[0:3]
            bi = (pos[i-1])[0:3] + "_" + uni
            tri = (pos[i-2])[0:3] + "_" + bi
                  
            if uni not in self.unilexicon:
                self.unilexicon[uni] = 0
            self.unilexicon[uni] += 1
            
            if bi not in self.bilexicon:
                self.bilexicon[bi] = 0
            self.bilexicon[bi] += 1
            
            if tri not in self.trilexicon:
                self.trilexicon[tri] = 0
            self.trilexicon[tri] += 1
            
            
            self.count += 1
            
    def finish(self):
        f = open("pos.txt","w")
        for word in sorted(self.unilexicon.keys()):
            f.write(word + "\t" + str(self.unilexicon[word]) + "\t" + str(self.unilexicon[word] / float(self.count)) + "\n")
        for word in sorted(self.bilexicon.keys()):
            f.write(word + "\t" + str(self.bilexicon[word]) + "\t" + str(self.bilexicon[word] / float(self.count)) + "\n")
        for word in sorted(self.trilexicon.keys()):
            f.write(word + "\t" + str(self.trilexicon[word]) + "\t" + str(self.trilexicon[word] / float(self.count)) + "\n")
        
            
if __name__ == "__main__":
    POSPP()
