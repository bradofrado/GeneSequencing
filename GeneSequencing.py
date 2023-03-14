#!/usr/bin/python3

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		self.init(seq1, seq2)

		for i in range(len(self.m)):
			for j in range(len(self.m[i])):
				insert = self.check_insert(i, j)
				delete = self.check_delete(i, j)
				match = self.check_match(i, j)
				minv = min(insert, delete, match)
				self.set_val(i, j, minv)
				self.set_pointer(i, j, self.get_pointer(i, j, insert, delete, match))
				

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = self.m[-1][-1];
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
	
	def init(self, seq1, seq2):
		self.m = [[None] * len(seq1) for i in range(len(seq2))]
		self.pointers = [[None] * len(seq1) for i in range(len(seq2))]
		self.seq1 = seq1
		self.seq2 = seq2


	def check_insert(self, i, j):
		if j - 1 < 0:
			return float('inf')
		
		return self.m[i][j-1] + INDEL

	def check_delete(self, i, j):
		if i - 1 < 0:
			return float('inf')
		
		return self.m[i - 1][j] + INDEL

	def check_match(self, i, j):
		last = self.m[i-1][j-1] if i-1 >= 0 and j-1 >= 0 else float('inf')
		
		val = SUB
		if self.seq1[i] == self.seq2[j]:
			if last == float('inf'):
				last = 0
			val = MATCH
		
		return val + last

	def set_val(self, i, j, val):
		self.m[i][j] = val
		pass
	def set_pointer(self, i, j, val):
		self.pointers[i][j] = val
		pass
	def get_pointer(self, i, j, insert, delete, match):
		if insert <= delete and insert <= match:
			return i,j-1
		if delete < insert and delete <= match:
			return i-1, j
		if match < delete and match < insert:
			return i-1, j-1
		Exception("ma")
