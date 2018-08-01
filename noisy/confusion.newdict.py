import os, sys, re, csv
import math, random
import pandas as pd
from collections import defaultdict

from scipy.stats import pearsonr, spearmanr

from pprint import pprint
from tqdm import tqdm


random.seed(123)
	
class word_type:

	def __init__(self, ortho, phones, freq, skeleton):
		self.ortho = ortho
		self.phones = phones
		self.freq = int(freq)
		self.skeleton = skeleton

	def __repr__(self):
		return self.ortho

	def __str__(self):
		return self.ortho

def read_confusion_matrix(f = 'c_initial.csv', f_dir = '../confusion/'):
	cm = {}
	df = pd.read_csv(f_dir + f, sep = '\t')
	stimulus = df.columns.values[0]
	resp = df.columns.values[1:]
	
	for i, x in enumerate(df[stimulus]):
		d = {}
		for j, y in enumerate(resp):
			d[y] = df[y][i] + 1

		sum_v = sum(d.values())
		
		d = dict((k, v / sum_v) for k, v in d.items())
		assert sum(d.values()) > .999 and sum(d.values()) < 1.001
		assert len(d) == len(resp)
		cm[x] = d
	return cm

def build_skeleton(word, Vs):
	# two passes
	skeleton = ['C' for _ in word]

	for i, p in enumerate(word):
		if p in Vs:
			skeleton[i] = 'V'
	# C for inital C
	# F for final c
	for i, p in enumerate(word):
		if i == 0 or skeleton[i] == 'V' or p == 'h':
			continue
		
		if p == 'G' or p == 'Z' or i == len(word) - 1:
			skeleton[i] = 'F'
		elif skeleton[i - 1] == 'V' and skeleton[i + 1] == 'C':
			skeleton[i] = 'F'
	
	for i, p in enumerate(word):
		if i == 0 or skeleton[i] != 'V':
			continue

		if i == len(word) - 1 or skeleton[i + 1] == 'C':
			skeleton[i] = 'W'

	skeleton = ''.join(skeleton)
	return skeleton

def read_nd(Vs, f = '../newdic.txt'):
	def syllabic_cons(m):
		return '^' + m.group(1).lower()
	t = {}
	with open(f) as rf:
		reader = csv.reader(rf, delimiter = '\t')
		for line in reader:
			phones, stress, main_stress, ortho, freq, pos = line
			phones = re.sub(r'x', '^', phones)
			phones = re.sub(r'\|', 'I', phones)
			phones = re.sub(r'R', 'X', phones)
			
			phones = re.sub(r'([MNL])', syllabic_cons, phones)
			skeleton = build_skeleton(phones, Vs)
			t[ortho] = word_type(ortho, phones, freq, skeleton)
	return t

def p_confuse(word_1, word_2, iCs, fCs, cVs, oVs):
	assert type(word_1) == word_type
	assert type(word_2) == word_type
	assert word_1.skeleton == word_2.skeleton
	skeleton = word_1.skeleton
	# 1 = stimulus
	# 2 = response
	total_p = 1
	for i, (c, phone_1, phone_2) in enumerate(zip(skeleton, word_1.phones, word_2.phones)):
		try:
			if c == 'V':
				p = cVs[phone_1][phone_2]
			elif c == 'W':
				p = oVs[phone_1][phone_2]
			elif c == 'F':
				p = fCs[phone_1][phone_2]
			elif c == 'C':
				p = iCs[phone_1][phone_2]
			total_p *= p
		except KeyError as k:
			print(c in cVs, c in oVs, c in fCs, c in iCs)
			print(i, skeleton, word_1.phones, word_2.phones)
			print(c, phone_1, phone_2, word_1, word_2)
			print(len(skeleton), len(word_1.phones))
			sys.exit()
		
	return total_p


def add_coca_freqs(nd, f = '../cf.txt'):
	cf = {}
	with open(f) as rf:
		reader = csv.reader(rf, delimiter = '\t')
		for line in tqdm(reader):
			word, count = line
			cf[word] = math.log(int(count))
			#cf[word] = int(count)
	
	for w in nd.copy():
		if w in cf:
			nd[w].freq = cf[w]
		else:
			del nd[w]
	return nd


def scramble_freqs(skeleton_to_words):
	length_to_freqs = defaultdict(list)
	for skeleton, words in skeleton_to_words.items():
		for w in words:
			length_to_freqs[len(skeleton)].append(w.freq)

	for l in length_to_freqs:
		random.shuffle(length_to_freqs[l])

	stw_copy = skeleton_to_words.copy()
	for skeleton, words in stw_copy.items():
		for w in words:
			w.freq = length_to_freqs[len(skeleton)].pop(0)

	return stw_copy

def log2(p):
	return math.log(p) / math.log(2)

def mutual_info(words, c_initial, c_final, v_closed, v_open, scramble = False, verbose = False):
	freq_sum = sum([w.freq for w in words])
	word_priors = dict((w.ortho, w.freq / freq_sum) for w in words)
	if scramble:
		words_as_list = list(word_priors.keys())
		priors_as_list = list(word_priors.values())
		random.shuffle(words_as_list)
		random.shuffle(priors_as_list)
		word_priors = dict((w, f) for w, f in zip(words_as_list, priors_as_list))

	H_y_x = 0
	marginals = {}
	for word_1 in words:
	
		# p(x)
		prior_p = word_priors[word_1.ortho]

		for word_2 in words:
			# p(y|x)
			confuse_p = p_confuse(word_1, word_2, c_initial, c_final, v_closed, v_open)
			H_y_x += prior_p * confuse_p * -log2(confuse_p)	
			# p(y)
			marginals[word_2.ortho] = marginals.get(word_2.ortho,0) + (prior_p * confuse_p)

	H_y = sum([p * -log2(p) for p in marginals.values()])
	I_X_Y = H_y - H_y_x
	if verbose:
		print('H(Y) =', H_y)
		print('H(Y|X) =', H_y_x)
		print('I(X;Y) =', I_X_Y)
	return I_X_Y

def output_data(skeleton_to_words, c_initial, c_final, v_closed, v_open, f = 'by_word.output.txt'):
	with open(f, 'w') as wf:
		writer = csv.writer(wf, delimiter = '\t')
		writer.writerow(['skeleton', 'length', 'word', 'p.as.self', 'freq'])
		for skeleton, words in skeleton_to_words.items():
			for w in words:
				confuse_p = p_confuse(w, w, c_initial, c_final, v_closed, v_open)
				writer.writerow([skeleton, len(skeleton), w.ortho, confuse_p, w.freq])

if __name__ == '__main__':
	l_type = 'dutch_'
	l_type = ''
	c_initial = read_confusion_matrix('{0}c_initial.csv'.format(l_type))
	c_final = read_confusion_matrix('{0}c_final.csv'.format(l_type))
	v_closed = read_confusion_matrix('{0}v_closed.csv'.format(l_type))
	v_open = read_confusion_matrix('{0}v_open.csv'.format(l_type))

	nd = read_nd(v_closed)
	nd = add_coca_freqs(nd)
	
	### organize words by their CV skeletons
	skeleton_to_words = defaultdict(set)
	for i, w in enumerate(sorted(nd.values(), key = lambda w : w.freq, reverse = True)):
		if len(w.skeleton) < 3:
			continue
		# get rid of skeletons that have consonant clusters
		if re.search('(CC|FF|[VW][VW]|[CF]{3})', w.skeleton):
			continue

		if len(skeleton_to_words[w.skeleton]) < 30:
			skeleton_to_words[w.skeleton].add(w)
	
	for skeleton, words in skeleton_to_words.copy().items():
		if len(words) < 5:
			del skeleton_to_words[skeleton]
			pass
	# H(Y|X) = Ex,y p(x,y) log p(y|x)
	#		= Ex,y p(y|x)P(x) log p(y|x)
	# H(Y) = Ey p(y) log p(y)
	
	
	output_data(skeleton_to_words, c_initial, c_final, v_closed, v_open, f =  '{0}by_word.output.txt'.format(l_type))
	sys.exit()
	
	Ps = []
	total_words = sum([len(words) for words in skeleton_to_words.values()])
	total_skeletons = len(skeleton_to_words) 
	print(total_words)
	print(total_skeletons)
	I_X_Ys = []
	for i in tqdm(range(1)):
		if i != -10:
			stw_scrambled = skeleton_to_words.copy()

		else:
			stw_scrambled = scramble_freqs(skeleton_to_words)	
		ixy = []
		for skeleton, words in stw_scrambled.items():
			if i == 0:
				mi = mutual_info(words, c_initial, c_final, v_closed, v_open, scramble = False)
			else:
				mi = mutual_info(words, c_initial, c_final, v_closed, v_open, scramble = True)

			ixy.append(mi)
		if i == 0:
			true_mi = sum(ixy)
			Ps.append([sum(ixy), 'actual'])
		else:
			Ps.append([sum(ixy), 'scrambled'])
			I_X_Ys.append(sum(ixy))
	
	print()
	mean_scr_mi = sum(I_X_Ys)/ len(I_X_Ys)
	sd = (sum([(mi - mean_scr_mi) ** 2 for mi in I_X_Ys]) / len(I_X_Ys)) ** .5
	print(mean_scr_mi, sd)
	print(mean_scr_mi + (2 * sd))
	print(true_mi)
	print(sum([true_mi >= mi for mi in I_X_Ys]) / len(I_X_Ys))
	import matplotlib.pyplot as plt
	plt.hist(I_X_Ys, bins = 10)
	plt.show()

	
	with open('{0}out.txt'.format(l_type), 'w') as wf:
		writer = csv.writer(wf, delimiter = '\t')
		writer.writerow(['mutual.info', 'lex.type'])
		writer.writerows(Ps)