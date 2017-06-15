# fmpl

**fmpl** is a Fortran code for structure learning of discrete graphical models.
**fmpl** takes as input a n x p data matrix and models the underlying p-dimensional probability distribution with a undirected graphical model
of (pairwise interacting) discrete random variables. The output of the code is a set of p * (p - 1) / 2 scores measuring the statistical interaction
between each pair of variables. More specifically, the scores are computed as corrected Frobenius norm scores
(see [Ekeberg et al.](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.87.012707)) from the coupling parameters obtained
by performing a (multinomial) logistic regression of each variable on the remaining variables
(see [Ravikumar et al.](http://projecteuclid.org/download/pdfview_1/euclid.aos/1268056617)).
We used **fmpl** to predict chromatin-related protein-protein interactions
(details in our paper ["Epigenomic co-localization and co-evolution reveal a key role for 5hmC as a communication hub in the chromatin network of ESCs"](http://www.sciencedirect.com/science/article/pii/S2211124716000280)) and contacts located at interfaces of eukaryotic protein heterodimers (["Conservation of coevolving protein interfaces bridges prokaryote--eukaryote homologies in the twilight zone"](http://www.pnas.org/content/113/52/15018.full)).

# LICENSE (BSD 3 clause)

Copyright (c) 2015, 2016, Simone Marsili
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
