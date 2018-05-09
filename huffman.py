import os, sys, re
import numpy as np
from pprint import pprint


_ab = 'abcdefghijklmnopqrstuvwxyz'

def is_num(x):
	return type(x) not in (list, tuple, np.ndarray)


class huff_2:
	def __init__(self, l):
		self.l = l
		self.t = self.binary_tree(range(len(self.l)))
		self.d = {}
		self.traverse_binary_tree(self.t)

	def rsum(self, l):
		if is_num(l):
			return self.l[l]
		else:
			elements = [self.rsum(little_l) for little_l in l]
			return sum(elements)

	def split_by_two(self, l):
		st = sorted(l, key = lambda x : self.rsum(x))
		node_1 = st.pop(0)
		node_2 = st.pop(0)
		new_node = [node_1, node_2]
		st.append(new_node)
		return st

	def binary_tree(self, t):
		while len(t) != 2:
			t = self.split_by_two(t)
		return t

	def traverse_binary_tree(self, t, s = ''):
		if is_num(t[0]):
			self.d[t[0]] = s + _ab[0]
		else:
			self.traverse_binary_tree(t[0], s + _ab[0])
		
		if len(t) == 2:
			if is_num(t[1]):
				self.d[t[1]] = s + _ab[1]
			else:
				self.traverse_binary_tree(t[1], s + _ab[1])

class huff_n:

	def __init__(self, p, n):
		self.l = p
		self.n = n
		self.t = self.build_tree(range(len(self.l)))
		self.d = {}
		self.traverse_tree(self.t)

	def rsum(self, l):
		if is_num(l):
			return self.l[l]
		else:
			elements = [self.rsum(little_l) for little_l in l]
			return sum(elements)

	def build_tree(self, t):
		while len(t) > self.n:
			t = sorted(t, key = lambda x : self.rsum(x))
			new_node = []
			for i in range(min(self.n, len(t))):
				new_node.append(t.pop(0))
			t.append(new_node)
		return t

	def traverse_tree(self, t, s = ''):
		for i in range(len(t)):
			if is_num(t[i]):
				self.d[t[i]] = s + _ab[i]
			else:
				self.traverse_tree(t[i], s + _ab[i])

np.random.seed(42)
p = np.random.randint(100, size = 10)
q = [1, 2, 3]


h = huff_2(q)
print(h.l)
print(h.t)
print(h.d)
h2 = huff_n(q, 3)
print(h2.t)
print(h2.d)
#print(h.d)
#pprint(sorted(h.d.values()))
