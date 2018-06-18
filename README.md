# bennu
Code and display software used for this paper on Bennu  
https://arxiv.org/abs/1804.02058

The mass-spring model code is built upon the N-body code rebound, 
see https://github.com/hannorein/rebound

The Bennu code example is in myexamples/bennu2/
and uses the Bennu shape model 
by Nolan, M.C., et al.  2013,
Shape model and surface properties of the OSIRIS-REx target Asteroid (101955)
Bennu from radar and light-curve observations.
Icarus 226, 629-640.

The bennu2 directory also contains parameter files used in the paper.

Associated jupyter notebooks for looking at outputs is in myexamples/pylab/  
(though I tend to move code outputs into subdirectories before running it).
The bennu_spec notebook lets you look at the spectrum of normal mode excitation.

Slightly edited rebound code is in src/

Spring routines (not part of rebound) are in src_spring/


