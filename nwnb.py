
# Attempt to make this function numba ready
import itertools
import scipy.special
import numba
import numpy as np
# A function for making a matrix of zeroes
#@numba.jit(nopython=True)
def nb_nw(seq1, seq2):
    """
    def zeros(rows, cols):
        # Define an empty list
        retval = [np.int64(x) for x in range(0)]
        # Set up the rows of the matrix
        for x in range(rows):
            # For each row, add an empty list
            retval.append([np.int64(x) for x in range(0)])
            # Set up the columns in each row
            for y in range(cols):
                # Add a zero to each column in each row
                retval[-1].append(0)
        # Return the matrix of zeros
        return retval"""

    # A function for determining the score between any two bases in alignment
    def match_score(alpha, beta, gap_penalty = -1,match_award = 1, mismatch_penalty = -1):
        if alpha == beta:
            return match_award
        elif alpha == '-' or beta == '-':
            return gap_penalty
        else:
            return mismatch_penalty

    def nb_needleman_wunsch(seq1, seq2, gap_penalty = -1,match_award = 1, mismatch_penalty = -1):

        # Store length of two sequences
        n = len(seq1)  
        m = len(seq2)
        
        # Generate matrix of zeros to store scores
        score = np.zeros((m+1, n+1))
       
        # Calculate score table
        
        # Fill out first column
        for i in range(0, m + 1):
            score[i][0] = gap_penalty * i
        
        # Fill out first row
        for j in range(0, n + 1):
            score[0][j] = gap_penalty * j
        
        # Fill out all other values in the score matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate the score by checking the top, left, and diagonal cells
                match = score[i - 1][j - 1] + match_score(seq1[j-1], seq2[i-1])
                delete = score[i - 1][j] + gap_penalty
                insert = score[i][j - 1] + gap_penalty
                # Record the maximum score from the three possible scores calculated above
                score[i][j] = max(match, delete, insert)
        
        # Traceback and compute the alignment 
        
        # Create variables to store alignment
        align1 = ""
        align2 = ""
        
        # Start from the bottom right cell in matrix
        i = m
        j = n
        
        # We'll use i and j to keep track of where we are in the matrix, just like above
        while i > 0 and j > 0: # end touching the top or the left edge
            score_current = score[i][j]
            score_diagonal = score[i-1][j-1]
            score_up = score[i][j-1]
            score_left = score[i-1][j]
            
            # Check to figure out which cell the current score was calculated from,
            # then update i and j to correspond to that cell.
            if score_current == score_diagonal + match_score(seq1[j-1], seq2[i-1]):
                align1 += seq1[j-1]
                align2 += seq2[i-1]
                i -= 1
                j -= 1
            elif score_current == score_up + gap_penalty:
                align1 += seq1[j-1]
                align2 += '-'
                j -= 1
            elif score_current == score_left + gap_penalty:
                align1 += '-'
                align2 += seq2[i-1]
                i -= 1

        # Finish tracing up to the top left cell
        while j > 0:
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        while i > 0:
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1
        
        # Since we traversed the score matrix from the bottom right, our two sequences will be reversed.
        # These two lines reverse the order of the characters in each sequence.
        align1 = align1[::-1]
        align2 = align2[::-1]
        
        return (align1, align2)
        
    a1, a2 = nb_needleman_wunsch(seq1, seq2)
    assert len(a1) == len(a2)
    tot = 0
    for i in range(len(a1)):
        if a1[i] != a2[i]:
            tot += 1

    return tot

@numba.jit(nopython=True, parallel=False)
def distance_vec(dvec, indices, seqs, nb_metric, *args):
    for veci in numba.prange(len(indices)):
        si = seqs[indices[veci, 0]]
        sj = seqs[indices[veci, 1]]
        d = nb_metric(si, sj, *args)
        dvec[veci] = d

def nb_pairwise_sq(seqs, nb_metric, *args):
    """Calculate distance between all pairs of seqs using metric
    and kwargs provided to nb_metric. Will use multiprocessing Pool
    if ncpus > 1.

    nb_metric must be a numba-compiled function

    Parameters
    ----------
    seqs : list
        List of sequences provided to metric in pairs.
    metric : numba-compiled function
        A distance function of the form
        func(seq1, seq2, **kwargs)
    **kwargs : keyword arguments
        Additional keyword arguments are supplied to the metric.

    Returns
    -------
    dvec : np.ndarray, length n*(n - 1) / 2
        Vector form of the pairwise distance matrix.
        Use scipy.distance.squareform to convert to a square matrix"""
    nb_seqs = numba.typed.List()
    for s in seqs:
        nb_seqs.append(s)

    dvec = np.zeros(int(scipy.special.comb(len(seqs), 2)))
    indices = np.zeros((int(scipy.special.comb(len(seqs), 2)), 2), dtype=np.int)
    for veci, ij in enumerate(itertools.combinations(range(len(seqs)), 2)):
        indices[veci, :] = ij

    distance_vec(dvec, indices, nb_seqs, nb_metric, *args)
    return dvec