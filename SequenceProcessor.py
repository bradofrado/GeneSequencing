import time
import math

from SequenceLoader import SequenceLoader
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
    
if __name__ == '__main__':
	loader = SequenceLoader()
	seqs = loader.load('test.txt')
	processor = SequenceProcessor(seqs)
	t, results, seqs = processor.process(True, 100)

	for result in results:
		for sub in result:
			if len(sub) == 0: 
				continue
			print('Assignment1: ' + sub['seqi_first100'] + ' Assignment2: ' + sub['seqj_first100'] + ' Score: ' + str(sub['align_cost']))
  
