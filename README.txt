*** High Performance Computing Project 2023 ***
Authors: Benjamin Gruber and Dominik Freinberger

Content:
- Overview of submitted files
- Description of the two main files 
    (.) 'allreduce.h'
    (.) 'benchmark.c'

---------------------
// FILES //
---------------------

-> README.txt
-> report.pdf
-> code
    -> 'allreduce.h'
    -> 'benchmark.c'
-> data
    -> pdf
        -> ex145
            -> 1x16.pdf
            -> 1x32.pdf
            -> 20x1.pdf
            -> 20x16.pdf
            -> 20x32.pdf
            -> 36x1.pdf
            -> 36x16.pdf
            -> 36x32.pdf
        -> ex2
            -> 1x16.pdf
            -> 1x32.pdf
            -> 20x1.pdf
            -> 20x16.pdf
            -> 20x32.pdf
            -> 36x1.pdf
            -> 36x16.pdf
        -> ex3
            -> 1x16.pdf
            -> 1x32.pdf
            -> 20x1.pdf
            -> 20x16.pdf
            -> 20x32.pdf
            -> 36x1.pdf
            -> 36x16.pdf
            -> 36x32.pdf
    -> csv
        -> ex145
            -> log_1x16_.out
            -> log_1x32_.out
            -> log_20x1_.out
            -> log_20x16_.out
            -> log_20x32_.out
            -> log_36x1_.out
            -> log_36x16_.out
            -> log_36x32_.out
        -> ex2
            -> log_1x16_.out
            -> log_1x32_.out
            -> log_20x1_.out
            -> log_20x16_.out
            -> log_20x32_.out
            -> log_36x1_.out
            -> log_36x16_.out
        -> ex3
            -> log_1x16_.out
            -> log_1x32_.out
            -> log_20x1_.out
            -> log_20x16_.out
            -> log_20x32_.out
            -> log_36x1_.out
            -> log_36x16_.out
            -> log_36x32_.out

Explanation:
------------
The 'report.pdf' file is the written report of all tests conducted.
The 'code' folder contains the code used to conduct the benchmarks.
The 'data' folder contains the raw data as collected on the Hydra system in two way,
either as pdf tables or as raw csv tables in the respective 'pdf' and 'csv' folders.
Due to time limits, the benchmarks for exercises 1, 4 and 5 were conducted in one run for 
each configuration whereas the tests for exercises 2 and 3 where conducted in two separate
runs.
For each of these three runs (1-4-5, 2 and 3) and for each NODES x PROCESSES configuration
the output of the algorithm is listed in a separate table. Note that for exercise 2, the 
configuration '36x32' did not finish in time and is therefore excluded.
Each table lists the following columns:
-ver:   the version of the algorithm (r...reference, 1...Exercise 1, ..., 5...Exercise 5) 
-size:  the total number of processes (as computed by nodes x tasks)
-count: the problem size (i.e. number of elements)
-b:     the block size
-avg:   the average runtime (depending on problem size, multiple runs where performed)
-std:   the standard deviation of the runtimes
        (for count = 10.000.000 the std is 0 as only one single run was performed each)
-med:   the median runtime
-min:   the minimum runtime
-lower: the lower bound of the 95% confidence interval
-upper: the upper bound of the 95% confidence interval

-------------------------
// DESCRIPTION OF CODE //
-------------------------

'allreduce.h'

Headerfile which contains all 5 custom Allreduce implementation according to the exercise sheet in
ascending order. The datatype to be used in the runs can be specified on top of the file, the default
is:
    #define TUW_TYPE MPI_FLOAT
    typedef float tuwtype_t;

'benchmark.c'

The program executes all algorithms contained in the allreduce.h file in one run for multiple 
problem sizes and block sizes. A correctness test is performed for the largest problem size
by comparing the result to the MPI_Allreduce implementation. A 'statistics' function is defined
which is executed each round by the MPI rank 0 process. The benchmark loops over various problem
sizes, block sizes and, depending on the problem size, does so multiple times to gather statistics.

Compilation was done with 

>> mpicc benchmark.c -o benchmark -O3 -lm
