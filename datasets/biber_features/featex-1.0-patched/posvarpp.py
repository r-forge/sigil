#!/usr/bin/env python

from preprocessor import Preprocessor
from input import *
import utils
from array import array

class POSPP(Preprocessor):
    
    def getName(self):
        return "POSVariancePreProcessor"
    
    def getAvailableOptions(self):
        return 'c:',['columns=']
    
    def parseOpts(self,opts):
        for (opt,arg) in opts:
            if opt == '-c' or opt == '--column':
                self.column = int(arg)
    
    def __init__(self):
        self.readNames()
        self.count = 0
        self.column = 1
        Preprocessor.__init__(self)
        
    
    def dis(self,pos):
        return pos[0:3]
    
    def readNames(self):
        f = open("posnames.txt")
        self.posnames = [l.strip() for l in f]
        self.posdict = {}
        for i in range(len(self.posnames)):
            self.posdict[self.posnames[i]] = i
        f.close()
    
    def run(self):
        self.counts = [array("d",[0.0 for j in range(len(self.files))]) for i in range(len(self.posnames))]
        self.filecount = 0
        for f in self.files:
            self.process(f)
    
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
            self.counts[self.posdict[uni]][self.filecount] += 1
            self.counts[self.posdict[bi]][self.filecount] += 1
            self.counts[self.posdict[tri]][self.filecount] += 1
            
            self.count += 1
        for x in self.posnames:
            self.counts[self.posdict[x]][self.filecount] /= float(len(pos)-3)
        self.filecount += 1
            
    def finish(self):
        f = open("pos_variances.txt","w")
        for i in range(len(self.posnames)):
            mean = utils.mean(self.counts[i])
            f.write(self.posnames[i] + "\t" + str(mean) + "\t" + str(utils.median(self.counts[i])) + "\t" + str(utils.variance(self.counts[i])) + "\t" + str(utils.moment(self.counts[i],mean,3)) + "\t" + str(utils.moment(self.counts[i],mean,4)) +  "\t" + str(len([x for x in self.counts[i] if x > 0])) + "\n")
            
if __name__ == "__main__":
    POSPP()
