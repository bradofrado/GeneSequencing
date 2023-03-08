import time
import math

from GeneSequencing import GeneSequencing

class SequenceProcessor:
    def __init__(self, seqs) -> None:
      self.seqs = seqs
      self.solver = GeneSequencing()
      pass
    def process(self, banded, alignLength):
      sequences = [ self.seqs[i][2] for i in sorted(self.seqs.keys()) ]
      processed_results = []
      start = time.time()

      for i in range(len(sequences)):
        jresults = []
        for j in range(len(sequences)):
          if(j < i):
            s = {}
          else:
            s = self.solver.align(sequences[i], sequences[j], banded=banded,
                            align_length=alignLength)
          jresults.append(s)
        processed_results.append(jresults)

      end = time.time()
      ns = end-start
      return ns, processed_results, self.seqs