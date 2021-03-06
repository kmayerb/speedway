"""
Some speed testing on pwseqdist 

"""
import sys
from fg_shared import _git
from os.path import join as opj
sys.path.append(opj(_git, 'pwseqdist'))
import pwseqdist as pw
from nw import *
from nwnb import nb_nw, nb_pairwise_sq
from seqs import mixed_seqs
import timeit
import numba
import numpy as np

# align some arbitrary seq
# Test that we can align all seqs without an error

def test_nested_for_loop_on_seqs(seqs, func):
    result = np.zeros((len(seqs), len(seqs)))
    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i >= j:
                result[i, j] = func(seqs[i], seqs[j])
    return result

if __name__ == "__main__":

    metrics = {'nw_parasail': pw.metrics.nw_hamming_metric,
               'nw_python' : py_nw,
               'nw_nocompile': nb_nw,
               'nw_numba': numba.jit(nb_nw, nopython=True)}

    metrics = {'nw_parasail': pw.metrics.nw_hamming_metric,
               'nw_numba': numba.jit(nb_nw, nopython=True)}

    parasail_res = pw.pairwise.apply_pairwise_sq(seqs=mixed_seqs, metric=pw.metrics.nw_hamming_metric, ncpus=1)
    
    n = 4
    for name in metrics:
        func = metrics[name]
        if not name == 'nw_numba':
            def test():
                res = pw.pairwise.apply_pairwise_sq(seqs=mixed_seqs*10, metric=func, ncpus=1)
            res = pw.pairwise.apply_pairwise_sq(seqs=mixed_seqs, metric=func, ncpus=1)
        else:
            def test():
                res = nb_pairwise_sq(mixed_seqs*10, func)
            res = nb_pairwise_sq(mixed_seqs, func)

        print("#####################################################")
        print("## %s ##" % name)
        x1 = timeit.timeit(test, number=n) / n
        print(f"\tTIME.IT METHOD: (AVE OF {n} RUNS: {round(x1, 3)} SECONDS")
        print(res[:10])
        print('%1.0f%%' % (100 * np.mean(res == parasail_res)))
