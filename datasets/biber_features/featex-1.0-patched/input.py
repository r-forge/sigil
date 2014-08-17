"""featex_v2 -- Linguistic Feature Extraction Framework in Python Reloaded.
File Input classes:

- InputReader: Handles basic file I/O, including compression.
"""

__version__ = "$Id: input.py 42 2007-06-30 10:20:10Z sm $"

import gzip
from os.path import basename, exists
import re

class InputReader:
    """Reads text from a file that can be gzip compressed.
    
    The name of the file should be given the constructor. The file is read
    when the read() method is called, and the text is then obtained by a call
    to getText(). In addition, a file id is extracted from the filename 
    (all text up to the first .), which is available though getID().
    
    If getText() is called before read(), this will trigger an implcit call to
    read(). The file is only read ONCE, all subsequent calls to read() will 
    have no effect.
    """
    def __init__(self,filename):
        self.__filename = filename
        self.__text = ''
        fn = basename(self.__filename)
        dotpos = fn.find('.')
        self.__id = fn[:dotpos]
        self.__isRead = False

    def getText(self):
        if not self.__isRead:
            self.read()
        return self.__text
    
    def getID(self):
        return self.__id
    
    def read(self):
        if not self.__isRead:
            if self.__filename[-2:] == 'gz':
                f = gzip.GzipFile(self.__filename) 
            else:
                f = open(self.__filename,'r')
            self.__text = f.read()
            f.close()
            self.__isRead = True
            
class CQPFormat:
    """Class that parses / extracts information from CQP-format files.
    
    CQP-format files are tab-separated files interleaved with SGML tags.
    A line in such a file contains either
        - nothing (will be ignored), or
        - a tab separated list of fields for one token, or
        - a opening or closing SGML tag (with optional attributes)
    """
    
    TAGREGEXP = re.compile("^</?([^> ]+)([^>]*)>$")
    
    def __init__(self,text):
        """Initialize the object by parsing the given text."""
        self.__text = text
        self.__stacks = {}
        self.__numColumns = 0
        self.__columns = []
        self.__pos = 0
        self.__annotations = {}
        
        self.__parse()
        
        
    def __parse(self):
        for line in self.__text.split('\n'):
            self.__parseLine(line)
            
    def __parseLine(self,line):
        stripped = line.strip()
        if stripped == '':
            return  # skip empty lines
        match = CQPFormat.TAGREGEXP.match(stripped)
        if match:
            self.__parseTag(stripped,match)
        else:
            self.__parseTab(line)
            
    def __parseTag(self,line,match):
        tag = match.group(1)
        if line[1] == "/":
            if tag not in self.__stacks:
                return
            if len(self.__stacks[tag]) != 0:
                start, attrs = self.__stacks[tag].pop()
                if tag not in self.__annotations:
                    self.__annotations[tag] = []
                self.__annotations[tag].append((start,self.__pos,attrs))
            else:
                raise ValueError("Closing a non-open tag " + tag)
    
        else:
            attrs = match.group(2).strip()
            if tag not in self.__stacks:
                self.__stacks[tag] = []
            self.__stacks[tag].append((self.__pos,attrs))
        
    
    def __parseTab(self,line):
        """Parses a line of tab-separated fields"""
        fields = line.split('\t')
        if self.__numColumns == 0:
            self.__numColumns = len(fields)
            for i in range(self.__numColumns):
                self.__columns.append([])
        
        elif self.__numColumns != len(fields):
            print line
            raise ValueError("Unequal number of columns!")
            
        for i in range(self.__numColumns):
            self.__columns[i].append(fields[i])
        
        self.__pos += 1
        
    def numColumns(self):
        return self.__numColumns
    
    def getColumn(self,columnId):
        return self.__columns[columnId]
    
    def getAnnotations(self,tag):
        return self.__annotations[tag]

class DepParseFormat:
    def __init__(self,text):
        self.__text = text
        self.__sentences = []
        self.__parse()
        print "processing ..."
        print self.__sentences[0]          

    def __parse(self):
        # Sentences are separated a blank line    
        sentences = self.__text.split('\n\n')
        for sentence in sentences[0:-1]:
            lines = sentence.split('\n')
            # The first line start with ## followed by the parsed sentence
            if lines[0][0:2] != '##':
                print sentence
                raise ValueError('Sentence does not start with comment!')
            deps = []
            for line in lines[1:]:
                deps.append(self.__parseLine(line))
            self.__sentences.append(deps)

    def __parseLine(self,line):
        fields = line.split('\t')
        return tuple(fields)
                   
