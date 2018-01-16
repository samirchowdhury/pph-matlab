# pph-matlab

Matlab code for computing path persistent homology of directed networks. 
The current implementation computes PPH up to homology dimension 1.

Chowdhury, S. and Mémoli, F., Persistent Path Homology of Directed Networks. To appear in SODA 2018.


### Example usage

A = rand(5).*(~eye(5)); % vertices appear at time 0, PPH cares about the entry of (directed) edges

pers = computePPH(A,max(max(A)))


### Further testing with cycle networks

load cycle_nets.mat 

A{3}

// A{3} should be a 5x5 matrix with entries 0 1 2 3 4 on top row
   
pers = computePPH(A{3},4); 

pers{2}

ans = 

   1  3

// More generally, A{n} contains a cycle network on (n+2) nodes, maximum weight = (n+1). 
By the characterization result in our paper, the 1-dimensional PPH of this object should be 
a single interval (1, ceil(n/2)). 
