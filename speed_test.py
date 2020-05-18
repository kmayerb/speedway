"""
Some speed testing on pwseqdist 

"""
import pwseqdist as pw
from nw import *
from nwnb import nb_nw
from seqs import mixed_seqs
import timeit

# align some arbitrary seq
# Test that we can align all seqs without an error
def test_nested_for_loop_on_seqs(seqs):
	result = list()
	for i in range(len(seqs)):
		for j in range(len(seqs)):
			if i >= j:
				x1,x2=needleman_wunsch(seqs[i], seqs[j])
				result.append((x1,x2))


if __name__ == "__main__":

	n = 3

	mixed_seqs = list(mixed_seqs + mixed_seqs + mixed_seqs)

	def wr1():
		test_nested_for_loop_on_seqs(seqs = mixed_seqs)

	print("#####################################################")
	print("## --- EDUCATIONAL PURE PYTHON NEEDLEMAN WUNSCH -- ##")
	x1 = timeit.timeit(wr1, number = n)/n
	print(f"METHOD:: test_nested_for_loop_on_seqs")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x1,2)} SECONDS")

	def wr3():
		pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = needleman_wunsch, ncpus = 1)

	x3 = timeit.timeit(wr3, number = n)/n
	print(f"METHOD: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = needleman_wunsch, ncpus = 1) (1 CPU)")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x3,2)} SECONDS")

	def wr4():
		pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = needleman_wunsch, ncpus = 6)
	
	x4 = timeit.timeit(wr4, number = n)/n
	print(f"METHOD: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = needleman_wunsch, ncpus = 8) (8 CPU)")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x4,2)} SECONDS")



	print("######################################")
	print("## --- PARASAIL NEEDLMAN WUNSCH -- ##")

	def wr5():
		pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_metric, ncpus = 1)
	
	x5 = timeit.timeit(wr5, number = n)/n
	print(f"METHOD:: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_metric, ncpus = 1)")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x5,2)} SECONDS")

	def wr6():
		pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_metric, ncpus = 6)

	x6 = timeit.timeit(wr6, number = n)/n
	print(f"METHOD:: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_metric, ncpus = 6")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x6,2)} SECONDS")


	print("######################################")
	print("## --- Parasail nw_hamming_metric- ##")

	def wr7():
		pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_hamming_metric, ncpus = 1)
	
	x7 = timeit.timeit(wr7, number = n)/n
	print(f"METHOD:: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_hamming_metric-, ncpus = 1)")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x7,2)} SECONDS")

	def wr8():
		pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_hamming_metric, ncpus = 6)

	x8 = timeit.timeit(wr8, number = n)/n
	print(f"METHOD: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_hamming_metric-, ncpus = 6)")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x8,2)} SECONDS")



	print("######################################")
	print("####.     NUMBA       ################")
	print("######################################")
	print("######################################")

	# DO NUMBA ON PREXISTING HAMMING (TRUNCATED SEQS TO BE OF SAME LENGTH)
	def wr0():
		mixed_seqs_5 = [s[0:5] for s in mixed_seqs]
		pw.numba_tools.nb_pairwise_sq(seqs = mixed_seqs_5, nb_metric = pw.numba_tools.nb_hamming_distance)

	x0 = timeit.timeit(wr0, number = n)/n
	print(f"METHOD: pw.numba_tools.nb_pairwise_sq(seqs = mixed_seqs_5, nb_metric = pw.numba_tools.nb_hamming_distance)")
	print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x0,2)} SECONDS")


	# 
	# def wr00():
	# 	mixed_seqs_5 = [s[0:5] for s in mixed_seqs]
	# 	pw.numba_tools.nb_pairwise_sq(seqs = mixed_seqs_5, nb_metric = nb_nw)

	# x00 = timeit.timeit(wr00, number = n)/n
	# print(f"METHOD: pw.numba_tools.nb_pairwise_sq(seqs = mixed_seqs_5, nb_metric = nb_needleman_wunsch)")
	# print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x00,2)} SECONDS")



