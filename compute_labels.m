function label_T = compute_labels(A,T)
% This function takes a cell array and adds labels of elementary paths
% Inputs:
% A is the network
% T is the initialized array

num_nodes = size(A,1);

  % Linear arrays are stored as rows of T. Label each linear array individually
  for ii=1:size(T,1)
    T{ii} = compute_individual_labels(A,num_nodes,T{ii},ii);
    end
  % Return T, each linear array labeled 
  label_T = T;  
  end
  

  
function label_individual_T = compute_individual_labels(A,num_nodes,Ti,k)
  %% This function labels a linear array.
  % Inputs:
  % A is the network
  % num_nodes is the number of nodes
  % Ti is the individual linear array
  % k is the parameter in n-permute-k. Should have k=1,2,3.
  
  % Initialize output
  label_individual_T = [];
  
  v = npermutek([1:num_nodes],k);
  %  remove consecutively repeated entries
  jj = 1;
  while (jj < k)
    [~,rows_to_kill] = pph_regularize(v(:,[jj,jj+1]));
    v = v(~rows_to_kill, :);
    jj = jj+1;
  end
    %    test that correct number of labels have been computed
    if (length(Ti) ~= size(v,1))
      disp('something wrong')
    end
      %      return labeled T
      v = num2cell(v,2);
      label_individual_T = v';
      end

  