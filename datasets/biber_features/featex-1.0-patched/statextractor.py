#!/usr/bin/env python
from extractor import Extractor
from input import InputReader, CQPFormat
import utils
import re

def compileRE(str):
        return re.compile(str, re.VERBOSE | re.IGNORECASE)
    

class StatExtractor(Extractor):
    """Extractor for some general statistics and some
    simple part of speech counts -- no real linguistic knowledge involved.
    """
    # GENERAL PATTERNS
    #######################################################
    WORDPOS = r"((?<=\s)(%s)/(%s)(?=\s))"
    POS = WORDPOS % ("\S+","%s")
    WORD = WORDPOS % ("%s","\S{3}")
    
    # two words
    WORDS_2 = WORD + "\s" + WORD
    
    #three words
    WORDS_3 = WORDS_2 + "\s" + WORD
    
    WORDS_4 = WORDS_3 + "\s" + WORD
    
    # HELPER PATTERNS, e.g. simplified tags such as N 
    #######################################################
    # forms of "have", including ones with negative particle attached
    H_HAVE = (POS % "VH\w") + "(\s" + (POS % "XX0") + ")?" 
    
    # auxiliary verbs
    H_AUX = (POS % "V[^V].")
    
    # past tense or past participled verb
    H_VBN = (POS % "V[BDHV][DN]")
    
    # any form of "be"
    H_BE = (POS % "VB.")
    
    # any kind of punctuation or start/end of sentence
    H_PUN = "(</?s>|" + (POS % "PU.") + ")"
    
    # any word (but not a sentence boundary!)
    H_ANY = POS % "\S{3}" 
    
    
    # List of all features
    DIRECT_FEATS = [ ]
    
    ## Pre-compile REs to speed things up a little
    DIRECT_FEATS = [(name, compileRE(regex)) for (name,regex) in DIRECT_FEATS]
    
    # features that need to be counted for calculations but are not directly
    # used.
    CALC_FEATS = [] 
   
    ## Pre-compile REs to speed things up a little
    CALC_FEATS = [(name, compileRE(regex)) for (name,regex) in CALC_FEATS]
   
             
    # the names of the features that are calculated on the basis of other
    # features 
    CALCULATED_FEATS = []
    
    # features to be counted on lemmatized from
    # list of tuples of the form (featName,wordset) 
    LEMMA_FEATS = []
    
    
    ### further possible non-Biber features
    # - all past participles (V.N)
    # - total subordinators (CJS) 
    
        
    def getName(self):
        return "StatExtractor"
    
    def getFeatureNames(self):
        return self.__featureNames
    
    def getAvailableOptions(self):
        short = ""
        long = []
        return short, long
        
    def __init__(self):
        self.__featureNames = [name for (name,re) in self.DIRECT_FEATS] + self.CALCULATED_FEATS
        Extractor.__init__(self)
        
    def disambiguatePOS(self,poslist):
        # For ambigious tags, use the more probable analysis (the first tag)
        return [x.split('-')[0] for x in poslist]
    
    def getWordsWithPOS(self,wordlist,poslist):
        if len(wordlist) != len(poslist):
            raise Error("Word- and POS-list are not of the same length!")
        return [word + "/" + pos for (word,pos) in zip(wordlist,poslist)]
    
    def extractWithREs(self,features,text):
        featDict = {}
        for (featName, featRE) in features:
            #print featName
            #print featRE.pattern
            matches = featRE.findall(text)
            #matches = re.findall("(" + featRE + ")", text, re.VERBOSE | re.IGNORECASE)
            #self.showMatches(featRE,text)
            count = len(matches)
            featDict[featName] = count
        return featDict
    
    def extractFromLemmatatizedForms(self,features,lemmaList):
        featDict = {}
        for (featName, words) in features:
            count = sum([1 for lemma in lemmaList if lemma in words])
            featDict[featName] = count
        return featDict
    
    def extractStatistics(self,cqpf):
        featDict = {}
        wordList = cqpf.getColumn(0)
        posList = cqpf.getColumn(1)
        lemmaList = cqpf.getColumn(2) # use lemmatized forms as tokens
        
        sentenceLengths = [end-start for (start,end,arg) in cqpf.getAnnotations("s")]
        featDict.update(utils.getStats("sent_length", sentenceLengths))
        
        wordLengths = []
        for i in range(len(wordList)):
            if posList[i] not in ["PUN", "PUL", "PUR", "PUQ"]:
                wordLengths.append(len(wordList[i]))
        featDict.update(utils.getStats("word_length", wordLengths))
        longWordLengths = [x for x in wordLengths if x > 8]
        featDict["long_words"] = len(longWordLengths) / float(len(wordLengths))
        featDict.update(utils.getStats("long_word_length", longWordLengths))
        
        
        
        
        tokens = len(lemmaList)
        types = len(frozenset(lemmaList))
        
        return featDict
    
    
    
    # calculate some features after extraction
    # the passed dictionary will be changed!
    def calculateFeats(self, featDict):
        pass
    
    def normalizeByLength(self,featDict,length):
        for feat in featDict:
            featDict[feat] = featDict[feat] / float(length) * 1000 
    
    def showMatches(self,pattern,text,num=10, context = 40):
        count = 0
        for match in pattern.finditer(text):
            start = match.start()
            end = match.end()
            print text[start-context:start] + " [ " + text[start:end] + " ] " + text[end:end+context]
            count += 1
            if count == num:
                return
    
    def process(self,file):
        feats = {}
        Extractor.process(self,file)
        ir = InputReader(file)
        ir.read()
        cqpf = CQPFormat(ir.getText())
        #words = ' '.join(cqpf.getColumn(0))
        #pos = ' '.join(self.disambiguatePOS(cqpf.getColumn(1)))
        lemma = cqpf.getColumn(2)
        sentences = cqpf.getAnnotations("s")
        wordpostmp = []
        for (start,end,attr) in sentences:
            wordpostmp.append('<s>')
            wordpostmp.extend(self.getWordsWithPOS(
                                cqpf.getColumn(0)[start:end],
                                self.disambiguatePOS(cqpf.getColumn(1)[start:end])))
            wordpostmp.append('</s> ')
        wordpos = ' '.join(wordpostmp)
        feats.update(self.extractWithREs(self.DIRECT_FEATS,wordpos))
        feats.update(self.extractWithREs(self.CALC_FEATS,wordpos))
        feats.update(self.extractFromLemmatatizedForms(self.LEMMA_FEATS,lemma))
        self.calculateFeats(feats)
        self.normalizeByLength(feats, len(lemma))
        feats.update(self.extractStatistics(cqpf))
        print feats
        return ir.getID(),feats
        
if __name__ == "__main__":
    StatExtractor()