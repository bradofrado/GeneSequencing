#!/usr/bin/python3

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class SeqNode:
	def __init__(self, val):
		self.val = val
		self.op = 0
		self.pointer = None
	def __str__(self):
		return self.op + ": " + str(self.val)
	
	def __repr__(self):
		return self.__str__()

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
				if i == 0 and j == 0:
					self.set_val(i, j, 0)
					self.set_op(i, j, MATCH)
					self.set_pointer(i, j, None)
					continue
				insert = self.check_insert(i, j)
				delete = self.check_delete(i, j)
				match, op = self.check_match(i, j)
				minv = min(insert, delete, match)
				self.set_val(i, j, minv)
				self.set_op(i, j, self.get_op(insert, delete, match, op))
				self.set_pointer(i, j, self.get_pointer(i, j, insert, delete, match))
				

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = self.m[-1][-1].val
		mod1, mod2 = self.get_modified_strings(seq1, seq2)
		alignment1 = '{}  DEBUG:({} chars,align_len={}{})'.format(mod1, 
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = '{}  DEBUG:({} chars,align_len={}{})'.format(mod2, 
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
	
	def init(self, seq1, seq2):
		self.m = [[None] * (len(seq1) + 1) for i in range(len(seq2) + 1)]
		self.pointers = [[None] * (len(seq1) + 1) for i in range(len(seq2) + 1)]
		self.seq1 = seq1
		self.seq2 = seq2

	def get_modified_strings(self, seq1, seq2):
		mod1 = '' + seq1
		mod2 = '' + seq2
		next = self.m[-1][-1]
		i = 1
		j = 1
		while next != None and i <= 100 and j <= 100:
			x = len(seq2) - i
			y = len(seq1) - j
			if next.op == 'insert':
				mod1 = insert(mod1, '-', y)
				i += 1
			elif next.op == 'delete':
				mod2 = replace(mod2, '-', x)
				j += 1
			#elif next.op == 'sub':
				#mod1 = replace(mod1, mod2[x], y)
			next = next.pointer
			i += 1
			j += 1

		return mod1, mod2	
				


	def check_insert(self, i, j):
		if j - 1 < 0:
			return float('inf')
		
		return self.m[i][j-1].val + INDEL

	def check_delete(self, i, j):
		if i - 1 < 0:
			return float('inf')
		
		return self.m[i - 1][j].val + INDEL

	def check_match(self, i, j):
		if i - 1 < 0 or j - 1 < 0:
			return float('inf'), 0
		last = self.m[i-1][j-1].val
		
		val = SUB
		if self.seq1[j - 1] == self.seq2[i - 1]:
			# if last == float('inf'):
			# 	last = 0
			val = MATCH
		
		return val + last, val

	def set_val(self, i, j, val):
		self.m[i][j] = SeqNode(val)
		pass
	def set_pointer(self, i, j, val):
		self.m[i][j].pointer = val
		pass
	def set_op(self, i, j, op):
		self.m[i][j].op = op

	def get_op(self, insert, delete, match, op = MATCH):
		if insert <= delete and insert <= match:
			return 'insert'
		if delete < insert and delete <= match:
			return 'delete'
		if match < delete and match < insert:
			return 'match' if op == MATCH else 'sub'
		Exception("ma")
	
	def get_pointer(self, i, j, insert, delete, match):
		if insert <= delete and insert <= match:
			return self.m[i][j-1]
		if delete < insert and delete <= match:
			return self.m[i-1][j]
		if match < delete and match < insert:
			return self.m[i-1][j-1]
		Exception("ma")

def replace(s, char, index):
	return s[:index] + char + s[index+1:]

def insert(s, char, index):
	if index >= len(s):
		return s + char
	return s[:index] + char + s[index:]