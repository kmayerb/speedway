# speedway
Speed Tests

## pwsd

```bash
python speed_test.py
```

```
#####################################################
## --- EDUCATIONAL PURE PYTHON NEEDLEMAN WUNSCH -- ##
METHOD:: test_nested_for_loop_on_seqs
	TIME.IT METHOD: (AVE OF 3 RUNS: 2.79 SECONDS
METHOD: pw.pairwise.apply_pairwise_sq (1 CPU)
	TIME.IT METHOD: (AVE OF 3 RUNS: 2.72 SECONDS
METHOD: pw.pairwise.apply_pairwise_sq (8 CPU)
	TIME.IT METHOD: (AVE OF 3 RUNS: 0.56 SECONDS
######################################
## --- PARASAIL NEEDLMAN WUNSCH -- ##
METHOD:: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_metric, ncpus = 1)
	TIME.IT METHOD: (AVE OF 3 RUNS: 0.29 SECONDS
METHOD:: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_metric, ncpus = 6
	TIME.IT METHOD: (AVE OF 3 RUNS: 0.12 SECONDS
######################################
## --- Parasail nw_hamming_metric- ##
METHOD:: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_hamming_metric-, ncpus = 1)
	TIME.IT METHOD: (AVE OF 3 RUNS: 0.29 SECONDS
METHOD: pw.pairwise.apply_pairwise_sq(seqs = mixed_seqs, metric = pw.metrics.nw_hamming_metric-, ncpus = 6)
	TIME.IT METHOD: (AVE OF 3 RUNS: 0.12 SECONDS
######################################
####.     NUMBA       ################
######################################
######################################
METHOD: pw.numba_tools.nb_pairwise_sq(seqs = mixed_seqs
```