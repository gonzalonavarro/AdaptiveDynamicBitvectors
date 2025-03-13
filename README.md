# AdaptiveDynamicBitvectors

This is a simple C library implementing the techniques in the paper

Gonzalo Navarro.
Adaptive Dynamic Bitvectors.
Proc. SPIRE'24, pages 204-217

which we hope you cite if using it for your research.

You do not need any auxiliary library to install; just run "make".

Apart from the bitvectors described in the paper, which are handled in the
files containing "BV" in their name, analogous code is provided for handling
dynamic arrays of integers, where you can just insert, delete, and access.
Those are handled in the files containing "Id" in their name.

Apart from the basic operations, the bitvectors implement hybridNext and
hybridNext0, which look for the next 1/0 from a given position. While this
can be done by composing rank and select, the given implementations are much
faster.

Parameterization
----------------

The most important parameter is Theta, for which I recommend the following
depending on the values of n and q:

If 1/q >= 0.1, use theta = 0.1.
Else, if 1/q <= 0.0001, use theta = 0.01.
Else, use theta = 0.01 if n <= 2^22 , and theta = 0.001 otherwise.

If you don't know 1/q even approximately, theta = 0.01 is relatively safe.

To change b or gamma, modify MaxBlockWords or Gamma in leafBV.c/leafId.c

To change theta/alpha/epsilon, modify Theta/Alpha/Epsilon in hybridBV.c/hybridId.

Other inner parameters can also be modified in hybridBV.c/hybridId.c


Examples of use
---------------

Files access.c, rank.c, select.c, and memory.c are the basic building blocks
of all the experiments, mixing those queries with updates for the desired
values of n and q, and measuring times. Execute without parameters to see their
usage. You can also build on them as examples on how to use the operations.

The file main.c performs more basic tests on the operations.
