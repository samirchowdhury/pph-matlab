# pph-matlab

Matlab code for computing path persistent homology of directed networks. 
The current implementation computes PPH up to homology dimension 1.

Chowdhury, S. and Mémoli, F., Persistent Path Homology of Directed Networks. To appear in SODA 2018.


### Example usage

A = rand(5).*(~eye(5)); % vertices appear at time 0, PPH cares about the entry of (directed) edges

pers = computePPH(A,max(max(A)))
