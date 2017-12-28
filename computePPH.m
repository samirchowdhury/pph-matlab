function pers = computePPH(A,t)
% This function computes 1-dimensional persistent path homology of a network
% Inputs:
% A is a square matrix 
% t is a threshold used to add all paths at some time  
% dim is the dimension we compute up to. default = 1
  dim = 1;
  
  %% Initialize linear arrays
  %  T is an empty collection of linear arrays 
  T = cell(dim+2,1);
  % nn = number of nodes
  nn = size(A,1);
  %  Compute number of elementary paths in each dimension. Then create linear arrays 
  %  with that many columns. Want to have: 
  %  T_labels for labels 
  %  T_allow_times for allow times
  %  T_entry_times for entry times.

  for ii = 1:dim+2
    ep_size = nn*((nn-1)^(ii-1));  
    T{ii} = zeros(1,ep_size);
    end
   
  % Pass A and T to a function to input labels
  T_labels = compute_labels(A,T);

  % Pass A and T_labels to a function to input allow times
  [T_labels, T_allow_times] = compute_times(A,T_labels,t);

  %  Compute barcode
  pers = compute_pers(A,T_labels,T_allow_times);

  end
  