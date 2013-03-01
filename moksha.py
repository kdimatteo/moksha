from Bio import SeqIO
from fuzzywuzzy import fuzz
import os
import re
import subprocess
import sys
import time
from compiler.ast import Break


'''
Rule #1: C = (C|T) 

Rule #2: Matches have a tolerance of up to 3, single digit mismatches
between pattern and dataset, not including C/T swaps (any character mismatches, not within the [C|T] set

Rule #3: Matches have a tolerance of up to 6 consecutive, single digit
mismatches between pattern and dataseterror.
'''


class Moksha:
    
    def __init__(self, pattern, datafile):
        self.input_pattern = pattern
        #input_pattern = 'TCCAAATGAAGTCATTATCAAA'
        self.data_file = datafile
        #datafile = 'opuntia.fasta'
        handle = open(self.data_file, "rU")
        
        for record in SeqIO.parse(handle, "fasta"):
            self.brute(record)
        handle.close()
        print "----------------------------\nMOKSHA Complete"
        
    
    '''    
    def rule_two(self):   
        # /Users/kdimatteo/bin/agrep-2.04/agrep -3 'T[C|T][C|T]AAATGAAGT[C|T]ATTAT[C|T]AAA' ./data/*.txt 
        pattern = self.input_pattern.replace("C", "[C|T]")
        
        handle = open(self.data_file, "rU")
        for record in SeqIO.parse(handle, "fasta"):
            # agrep only works with files
            tmpfile = '/tmp/moksha_%s.txt' % time.time()
            h = open(tmpfile, 'w') #because you are too lazy to escape the string properly
            h.write(str(record.seq))
            h.close()
            
            cmd = "/Users/kdimatteo/bin/agrep-2.04/agrep -3 '%s' %s" % (pattern, tmpfile)
            print cmd
            #print subprocess.call(cmd, shell=True)  
           
            #p1 = subprocess.Popen([cmd], stdout=subprocess.PIPE)
            #print p1.communicate()[0]
            #p1.close()
    '''
    
        
        
    def brute(self, sequence_obj):
        print "processing %s" % sequence_obj.name
        # n can be q
        #s = str(sequence_obj.seq)
        #s = "somebaqjotimes"
        #pattern = "ban"
        
        s = str(sequence_obj.seq).strip()
        p = self.input_pattern
        
        #print "testing: %s" % sequence_obj.id
        
        L = [] #pattern for verification at end of the run 
        Q = [] #string from the sequence for verification

        digit_errors = 0 #mismatches, up to 3
        space_errors = 0 #spaced between errors, up to 6
        subs = 0
        seqIndex = 1
        
        i = 0 
        j = 0
        
        while i in range(0, len(s)):
            #print s[i]
            
            if len(Q) == len(p):
                print "==============================================\n"
                break
            
            #print "s: %s : %s  ?= p: %s : %s (e: %s) %s (%s)" % (i, s[i], j, p[j], digit_errors, "".join(L), len(L))
            if str(s[i]) == str(p[j]):
                L.append(p[j])
                Q.append(s[i])
                i = i + 1
                j = j + 1
                
            elif s[i] == 'T' and p[j] == 'C':
                L.append(p[j])
                Q.append(s[i])
                i = i + 1
                j = j + 1
                subs = subs + 1
            else:
                # no match
                if len(L) > 0:
                    if digit_errors <= 3:
                        digit_errors = digit_errors + 1
                        space_errors = space_errors + 1
                        L.append(p[j])
                        Q.append(s[i])
                        i = i + 1
                        j = j + 1
                        
                    else:
                        i = seqIndex
                        seqIndex = seqIndex + 1
                        j = 0
                        digit_errors = 0;
                        L = []
                        Q = []
                else:
                    i = i + 1
                    j = 0
                    # nothin in L anway
                    pass
                    
            #print i
                            
                    
            '''
            if space_errors >= 6:
                space_errors = 0
                digit_errors = 0
            '''    
                 
                    
        
        '''
        for i in s:
            fred = False
            
            j=0
            for j in range(0, len(pattern)):
                
              # j get added to test for success
              #  print "seq:%s :: pat:%s" % (i, j)
                if str(pattern[j]) == str(i):
                    L.append(pattern[j])
                    j = j + 1
                    break
                elif i == 'C' and pattern[j] == 'T':
                    L.append("Q")
                    j = j + 1
                    subs = subs + 1
                    print i, pattern[j]
        '''
        '''
        elif fred == False and i == 'C' and j == 'T':
            subs = subs + 1
            L.append('Q')
            break
        
        else :
            continue
         '''
            
        '''
        if errors <= 160:
            if j == i:
                L.append(j)
                match = True
            else:
                errors = errors +1 
                
                match = False
          
            if j == i:
               L.append(j)
            elif j == 'T' and i == 'C':
                L.append(j)
                subs = subs + 1
            else:
                L.append(j)
                errors += 1
          '''
                    
        #print "s>" , "".join(Q)
        #print "p>" , "".join(L)
        
        #print "<<" , self.input_pattern
        
        
        if "".join(L) == self.input_pattern:
            print "[@][SUCCESS] %s" % sequence_obj.description
            print "Found this in the sequence:", "".join(Q)
        else:
            pass
            #print "[!] Failure %s" % sequence_obj.description
        #print "subs: %s digit_errors : %s space_errors :%s " % (subs, digit_errors, space_errors)            
        #print "\n"
                
        ''''
        if i != pattern[j]:
            errors = errors + 1
            break
        else:
           print pattern[j]
           chunk.append(pattern[j])
           j = j + 1
        '''
            
            
    def agrep(self, sequence_obj):
        # agrep only works with files
        pattern = self.input_pattern.replace("C", "[C|T]")
        tmpfile = '/tmp/moksha_%s.txt' % time.time()
        h = open(tmpfile, 'w')
        h.write(str(sequence_obj.seq))
        h.close()
        cmd = "/Users/kdimatteo/bin/agrep-2.04/agrep -3 '%s' %s" % (pattern, tmpfile)
        p = os.popen(cmd)
        print p.read()
        if (len(p.read()) > 0):
            print sequence_obj.id

    def fuzzy_match(self, sequence_obj):
        # see http://seatgeek.com/blog/dev/fuzzywuzzy-fuzzy-string-matching-in-python
        m = fuzz.partial_ratio(self.input_pattern, str(sequence_obj.seq))
        if(m>90):
            print sequence_obj.id
        
                    
    def regex_match(self, sequence_obj):
        pattern = self.input_pattern.replace("C", "[C|T]")
        p = re.compile(pattern)
        match_strings = p.findall(str(sequence_obj.seq))
        if len(match_strings)>0:
            print "Found: %s in %s" % (", ".join(match_strings), sequence_obj.id)
            
            m = p.search(str(sequence_obj.seq))
            if m != None:   
                print "Found %s at [%s, %s] " % (m.group(), m.start(), m.end())
        

    
    
# python moksha.py TCCAAATGAAGTCATTATCAAA data/CCDS_nucleotide.current.fna    
def test():
    o = Moksha(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    print "begin Moksha"
    from timeit import Timer
    #o = Moksha(sys.argv[1], sys.argv[2])
    t = Timer("test()", "from __main__ import test")
    print t.timeit()




    