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
		self.d = 3
		self.MaxCharactersToAlign = align_length

		self.init(seq1, seq2)
		self.calc_edit_distance_banded() if self.banded else self.calc_edit_distance()
				

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		x = self.d if self.banded else -1
		if self.banded:
			while x > 0 and self.m[-1][x] == None:
				x -= 1
		if self.m[-1][x] == None:
			score = float('inf')
			alignment1 = 'No Alignment Possible';
			alignment2 = alignment1
		else: 
			score = self.m[-1][x].val
			mod1, mod2 = self.get_modified_strings(seq1, seq2)
			alignment1 = mod1
			alignment2 = mod2
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
	
	def init(self, seq1, seq2):
		self.seq1_length = min(len(seq1), self.MaxCharactersToAlign)
		self.seq2_length = min(len(seq2), self.MaxCharactersToAlign)
		len1 = 2*self.d if self.banded else self.seq1_length
		len2 = self.seq2_length
		self.m = [[None] * (len1 + 1) for i in range(len2 + 1)]
		self.pointers = [[None] * (len1 + 1) for i in range(len2 + 1)]
		self.seq1 = seq1
		self.seq2 = seq2

	def calc_edit_distance(self):
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
	
	def calc_edit_distance_banded(self):
		for i in range(len(self.m)):
			for j in range(self.d * 2 + 1):
				if j < self.d - i or j - self.d + i > len(self.seq1):
					continue
				elif i == 0 and j == self.d:
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


	def get_modified_strings(self, seq1, seq2):
		mod1 = '' + seq1[:self.seq1_length]
		mod2 = '' + seq2[:self.seq2_length]
		x = self.d if self.banded else -1
		if self.banded:
			while x > 0 and self.m[-1][x] == None:
				x -= 1
		next = self.m[-1][x]
		i = 0
		j = 0
		while next != None: #and i <= 100 and j <= 100:
			x = self.seq2_length - i
			y = self.seq1_length - j
			if next.op == 'insert':
				mod2 = insert(mod2, '-', x)
				j += 1
			elif next.op == 'delete':
				mod1 = insert(mod1, '-', y)
				i += 1
			#elif next.op == 'sub':
				#mod1 = replace(mod1, mod2[x], y)
			next = next.pointer
			i += 1
			j += 1

		return mod1, mod2	
				


	def check_insert(self, i, j):
		lim = self.d - i if self.banded else 0
		if j - 1 < lim or j - 1 < 0:
			return float('inf')
		
		return self.m[i][j-1].val + INDEL

	def check_delete(self, i, j):
		x = j + 1 if self.banded else j
		if i - 1 < 0 or x >= len(self.m[i]):
			return float('inf')
		
		return self.m[i - 1][x].val + INDEL

	def check_match(self, i, j):
		x = j if self.banded else j - 1
		lim = self.d - i + 1 if self.banded else 0
		if i - 1 < 0 or x < lim:
			return float('inf'), 0
		last = self.m[i-1][x].val
		
		val = SUB
		x = j - self.d + i if self.banded else j
		if self.seq1[x - 1] == self.seq2[i - 1]:
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
			x = j + 1 if self.banded else j
			return self.m[i-1][x]
		if match < delete and match < insert:
			x = j if self.banded else j - 1
			return self.m[i-1][x]
		Exception("ma")

def replace(s, char, index):
	return s[:index] + char + s[index+1:]

def insert(s, char, index):
	if index >= len(s):
		return s + char
	return s[:index] + char + s[index:]