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
		self.d = 3 if banded else 0
		self.MaxCharactersToAlign = align_length

		self.init(seq1, seq2)
		self.calc_edit_distance_banded() if banded else self.calc_edit_distance()
				

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
		self.seq1_length = len(seq1)
		if self.MaxCharactersToAlign < self.seq1_length:
			self.seq1_length = self.MaxCharactersToAlign
		self.seq2_length = len(seq2)
		if self.MaxCharactersToAlign < self.seq2_length:
			self.seq2_length = self.MaxCharactersToAlign

		self.len1 = 2*self.d + 1 if self.banded else self.seq1_length + 1
		self.len2 = self.seq2_length + 1
		self.m = [[None] * self.len1 for i in range(self.len2)]
		self.pointers = [[None] * self.len1 for i in range(self.len2)]
		self.seq1 = seq1
		self.seq2 = seq2

	def calc_edit_distance(self):
		self.m[0][0] = 0
		for i in range(1, self.len1):
			self.m[0][i] = i * INDEL
			self.pointers[0][i] = (0, -1, 0)
		for i in range(1, self.len2):
			self.m[i][0] = i * INDEL
			self.pointers[i][0] = (-1, 0, 1)
		for i in range(1, self.len2):
			for j in range(1, self.len1):
				diff = MATCH if self.seq1[j - 1] == self.seq2[i - 1] else SUB
				insert = self.m[i][j-1] + INDEL
				delete = self.m[i - 1][j] + INDEL
				match = self.m[i-1][j-1] + diff
				if insert <= delete and insert <= match:
					minv = insert
					pointer = (0, -1, 0)
				elif delete < insert and delete <= match:
					minv = delete
					pointer = (-1, 0, 1)
				else:
					minv = match
					pointer = (-1, -1, 2)
				self.m[i][j] = minv
				self.pointers[i][j] = pointer

	def calc_edit_distance_banded(self):
		startX = self.d
		seq1_len = len(self.seq1)
		cnt = 1
		self.m[0][startX] = 0
		for i in range(startX + 1, self.len1):
			self.m[0][i] = (i - startX) * INDEL
			self.pointers[0][i] = (0, -1, 0)
		endY = self.d + 1
		for i in range(1, endY):
			x = max(self.d - i, 0)
			self.m[i][x] = i * INDEL
			self.pointers[i][x] = (-1, 0 + cnt, 1)
		startJ = 0
		for i in range(1, self.len2):
			for j in range(startJ, self.len1):
				if i + j <= self.d or j - self.d + i > seq1_len:
					continue
				x = j - self.d + i
				diff = MATCH if self.seq1[x - 1] == self.seq2[i - 1] else SUB
				if j - 1 >= 0:
					insert = self.m[i][j-1] + INDEL
				else:
					insert = float('inf')
				if j + cnt < self.len1:
					delete = self.m[i - 1][j + cnt] + INDEL
				else:
					delete = float('inf')
				match = self.m[i-1][j-1 + cnt] + diff
				if insert <= delete and insert <= match:
					minv = insert
					pointer = (0, -1, 0)
				elif delete < insert and delete <= match:
					minv = delete
					pointer = (-1, 0 + cnt, 1)
				else:
					minv = match
					pointer = (-1, -1 + cnt, 2)
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

			op = next[2]
			if op == 0:
				mod2 = '-' + mod2
				mod1 = seq1[y] + mod1
				i -= 1
			elif op == 1:
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