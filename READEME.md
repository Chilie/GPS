# Solving Phase Retrieval via Graph Projection Splitting

This package is to accompany our arXiv paper for reproducible research. It is very primitive and some polishes are needed. 

As far as I know, the phase transition of GPS for Gaussian phase retrieval is very sharp. In my opinion, our projection-based method GPS is very insensitive to the initialization.  And this is the big difference between other nonconvex solvers. Without the help of good initialization, our algorithm outperforms other nonconvex solvers.  See the benefits and the numerical experiments in our accompany paper. 

This package also includes some test files to ease the usage of the package, and you can just download it and run the demo.m file in Matlab.  Other demo files include Fourier phase retrieval with TV minimization regularization. Combined with this kind of regularization, the quality of reconstruction improves.

Any suggests and comments are welcome!

Reference: [Solving Phase Retrieval via Graph Projection Splitting](https://arxiv.org/pdf/1910.08714.pdf)

Feel free to contact me at matliji@nus.edu.sg