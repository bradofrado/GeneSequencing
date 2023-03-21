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
		self.d = 3
		self.MaxCharactersToAlign = align_length

		self.init(seq1, seq2)
		self.calc_edit_distance()
				

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		x = self.d if self.banded else -1
		if self.banded:
			while x > 0 and self.m[-1][x] == None:
				x -= 1
		if self.m[-1][x] == None:
			score = float('inf')
			alignment1 = 'No Alignment Possible'
			alignment2 = alignment1
		else: 
			score = self.m[-1][x]
			mod1, mod2 = self.get_modified_strings(seq1, seq2)
			alignment1 = mod1[:100]
			alignment2 = mod2[:100]
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
	
	def init(self, seq1, seq2):
		self.seq1_length = min(len(seq1), self.MaxCharactersToAlign)
		self.seq2_length = min(len(seq2), self.MaxCharactersToAlign)
		self.len1 = 2*self.d + 1 if self.banded else self.seq1_length + 1
		self.len2 = self.seq2_length + 1
		self.m = [[None] * self.len1 for i in range(self.len2)]
		self.pointers = [[None] * self.len1 for i in range(self.len2)]
		self.seq1 = seq1
		self.seq2 = seq2

	def calc_edit_distance(self):
		startX = self.d if self.banded else 0
		startY = 0
		seq1_len = len(self.seq1)
		for i in range(self.len2):
			for j in range(self.len1):
				if self.banded and (j < self.d - i or j - self.d + i > seq1_len):
					continue
				if i == startY and j == startX:
					self.m[i][j] = 0
					continue
				insert = self.check_insert(i, j)
				delete = self.check_delete(i, j)
				match = self.check_match(i, j)
				if insert <= delete and insert <= match:
					minv = insert
					pointer = 0, -1
				elif delete < insert and delete <= match:
					minv = delete
					pointer = -1, (1 if self.banded else 0)
				else:
					minv = match
					pointer = -1, (0 if self.banded else -1)
				self.m[i][j] = minv
				self.pointers[i][j] = pointer

	def get_modified_strings(self, seq1, seq2):
		mod1 = ''
		mod2 = ''
		currY = len(self.m) - 1
		currX = self.d if self.banded else len(self.m[0]) - 1
		if self.banded:
			while currX > 0 and self.pointers[currY][currX] == None:
				currX -= 1
		next = self.pointers[currY][currX]
		i = 1
		j = 1
		while next != None:
			x = self.seq2_length - i
			y = self.seq1_length - j
			currX += next[1]
			currY += next[0]

			op = self.get_op(next)
			if op == 'insert':
				mod2 = '-' + mod2
				mod1 = seq1[y] + mod1
				i -= 1
			elif op == 'delete':
				mod1 = '-' + mod1
				mod2 = seq2[x] + mod2
				j -= 1
			else:
				mod1 = seq1[y] + mod1
				mod2 = seq2[x] + mod2
			next = self.pointers[currY][currX]
			i += 1
			j += 1

		return mod1, mod2	
				


	def check_insert(self, i, j):
		lim = self.d - i if self.banded else 0
		if j - 1 < lim or j - 1 < 0:
			return float('inf')
		
		return self.m[i][j-1] + INDEL

	def check_delete(self, i, j):
		x = j + 1 if self.banded else j
		if i - 1 < 0 or x >= len(self.m[i]):
			return float('inf')
		
		return self.m[i - 1][x] + INDEL

	def check_match(self, i, j):
		x = j if self.banded else j - 1
		lim = self.d - i + 1 if self.banded else 0
		if i - 1 < 0 or x < lim:
			return float('inf')
		last = self.m[i-1][x]
		
		val = SUB
		x = j - self.d + i if self.banded else j
		if self.seq1[x - 1] == self.seq2[i - 1]:
			val = MATCH
		
		return val + last
	
	def get_op(self, coords):
		x, y = coords
		if x == 0 and y == -1:
			return 'insert'
		if x == -1 and y == (1 if self.banded else 0):
			return 'delete'
		if x == -1 and y == (0 if self.banded else -1):
			return 'match' 