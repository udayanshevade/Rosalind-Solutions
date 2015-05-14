import timeit

#exercise 0
#import this

# WARM UP -------------------------------

#exercise 1 - basic arithmetic
def hypotenuse(a, b):
	return a**2 + b**2
#print (hypotenuse(839, 874))


#exercise 2 - string slicing
def string_slice(string, a, b, c, d):
	return string[a:b+1] + " " + string[c:d+1]

sample_string_1 = "aAw1CaTezu7AgbII6GbNRqgndJOg1ZYkNFoA53BDQDjMachetes\
CkgzOvzBsUxrC00w3aC3sDuXZn2nVDEn9nxL9UEka6iridh1nJBA\
Z9ghPyiKuCgZdoGlhR1IiQrgregaria82YT7n9XbsTxNl0NUho6i\
pW1hZD6xD5VkZasAcxzriGja57Ob1iesJku4NA."

#print (string_slice(sample_string_1, 43, 50, 126, 133))

#exercise 3 - practicing loops
def odd_sum(a, b):
	return sum([i for i in range(a,b+1) if i % 2 != 0])
#print (odd_sum(4172, 8175))

#exercise 4 - working with files
def read_file_stepped(path):
	text = open(path).readlines()
	for i, line in enumerate(text):
		if "\n" in line:
			line = line[:-1]
		if i % 2 == 1:
			print (line)

#sample_file = r"C:\Users\Udayan\Desktop\Rosalind\python_village\read_text_staggered.txt"
#(read_file_stepped(sample_file))


#exercise 5 - practicing with dictionaries
def word_freq(path):
	freq_graph = {}
	text = open(path).read()
	for word in text.split():
		if word not in freq_graph:
			freq_graph[word] = 1
		else:
			freq_graph[word] += 1
	for key in freq_graph:
		print (key, freq_graph[key])

#sample_path = r"C:\Users\Udayan\Desktop\Rosalind\python_village\dictionary_practice.txt"
#word_freq(sample_path)


# START ---------------------------------

#common fxns
def read_file(path):
	return open(path).read()

def GC_count(string):
	return (string.count("G") + string.count("C"))/len(string) * 100

#----------------------------------------

#exercise 1 - intro molecular biology

# Counts bases in a nucleic acid string
path1 = r"C:\Users\Udayan\Desktop\Rosalind\bio_stronghold\count_bases.txt"

def letter_freq(seq):
	nucs = ["A", "C", "G", "T"]
	for nuc in nucs:
		print (seq.count(nuc))

#letter_freq(read_file(path1))

#----------------------------------------

#exercise 2 - second nucleic acid

# Replaces every instance of a given base with another
path2 = r"C:\Users\Udayan\Desktop\Rosalind\bio_stronghold\replace_uracil.txt"

def replace_uracil(seq):
	return seq.replace("T", "U")

#print (replace_uracil(read_file(path2)))

#----------------------------------------

#exercise 3 - complementing a strand

# Returns the complementary strand of code to the given string
path3 = r"C:\Users\Udayan\Desktop\Rosalind\bio_stronghold\complement_strand.txt"

def strand_complement(seq):
	complements= {"A": "T",
				"G": "C",
				"C": "G",
				"T": "A"}
	nucleotides = []
	for base in seq:
		nucleotides.append(base)
	for i in range(0, len(nucleotides)):
		nucleotides[i] = complements[nucleotides[i]]
	complement = "".join(nucleotides)
	return complement[::-1]

#print (strand_complement(read_file(path3)))


#----------------------------------------

#exercise 4 - wascally wabbits

# Calculates the population of rabbits given a steady rate
def rabbit_season(n, k):
	rabbits = {1: 1, 2:1}
	for i in range(3, n+1):
		rabbits[i] = rabbits[i-1] + rabbits[i-2]*k
		del rabbits[i-2]
	return rabbits[n]

#print (rabbit_season(28, 4))

#----------------------------------------

path4 = r"C:\Users\Udayan\Desktop\Rosalind\bio_stronghold\gc_content.txt"

# Finds the percentage of a strand that is made of G and C bases
def measure_GC_content(file_path):
	s = open(file_path)
	line = s.readline().rstrip("\n")
	max_GC, max_id = 0, ""
	while line:
		current_string, current_id = "", line[1:]
		line = s.readline().rstrip("\n")
		while ">" not in line and line:
			current_string += line
			line = f.readline().rstrip("\n")
		if current_string and GC_count(current_string) > max_GC:
			max_GC, max_id = GC_count(current_string), current_id
	print (max_GC)
	print (max_id)

#(measure_GC_content(path4))

#----------------------------------------

# Finds the variation between different strands
path5 = r"C:\Users\Udayan\Desktop\Rosalind\bio_stronghold\pt_mutations.txt"
text = open(path5)
sample_string1 = text.readline().rstrip()
sample_string2 = text.readline().rstrip()

def hamming_distance(a,b):
	return sum([1 for (x,y) in zip(a,b) if x != y])
	# count = 0
	# for x,y in zip(a,b):
	# 	if x != y:
	# 		count += 1
	# return count

#print (hamming_distance(sample_string1, sample_string2))

#----------------------------------------

#----------------------------------------

# Translates a strand into its corresponding amino acid chain
path6 = r"C:\Users\Udayan\Desktop\Rosalind\bio_stronghold\translation.txt"
text = open(path6)
sample_string3 = text.readline().rstrip()

def translation(RNA_string):
	aa_table = {"UUU": "F", "UUC": "F",
				"UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
				"AUU": "I", "AUC": "I", "AUA": "I",
				"AUG" : "M",
				"GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
				"UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
				"CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
				"ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
				"GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
				"UAU": "Y", "UAC":"Y",
				"CAU": "H", "CAC": "H",
				"CAA": "Q", "CAG": "Q",
				"AAU": "N", "AAC": "N",
				"AAA": "K", "AAG": "K",
				"GAU": "D", "GAC": "D",
				"GAA": "E", "GAG": "E",
				"UGU": "C", "UGC": "C",
				"UGG": "W",
				"CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
				"AGU": "S", "AGC": "S",
				"GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
	protein = []
	for i in range(0, len(RNA_string)-1, 3):
		codon = RNA_string[i] + RNA_string[i+1] + RNA_string[i+2]
		if codon in ["UAA", "UAG", "UGA"]:
			return "".join(protein)
		protein.append(aa_table[codon])
	return "".join(protein)

print (translation(sample_string3))

#----------------------------------------

# To be continued

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

#----------------------------------------

