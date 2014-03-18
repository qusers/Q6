#!/bin/bash
# Remember to change here the number of nodes to what is available in the queue
# you're calling. Perhaps you'll also need to issue the --exclusive option.
sbatch -n 4 -t 02:00:00 -J testqdynp5 run_test_mpi.sh

