function [T_labels, T_allow_times] = compute_times(A,T_labels,t_end)
% Compute allow and entry times 
% Inputs:
% A is the network
% T_labels is a collection of linear arrays, one for each dimension. Each entry
% contains the label for an elementary path.
% t_end is a threshold used to add all paths at some time

% Initialize T_allow_times as T_labels. This gives the right dimensions. We will
% switch entries of T_allow_times one by one. 

T_allow_times = T_labels;



% Compute some useful constants
num_nodes = size(A,1);
num_dims = size(T_labels,1);
t_start = min(min(A));

% Throw error if t_start >= t_end
if (t_start >= t_end)
  error('start time greater than end time');
end
  

% 0-path convention: in the digraph filtration, we start with all nodes and add
% edges incrementally. So set 0-path entry times to be minimum of A. 
len = length(T_allow_times{1});
v = ones(1,len).*t_start;
T_allow_times{1} = num2cell(v);

% 1-paths are used for getting other paths. Compute their allow times and store
len = length(T_allow_times{2});
lab = T_labels{2};
  % Get labels into matrix format
  lab_vector = cell2mat(lab);
  lab_array = reshape(lab_vector,2,len)';
  %  convert labels into linear indices of A
  lab_indices = sub2ind(size(A), lab_array(:,1), lab_array(:,2));
    % Get 1-path allow times
    v = ones(1,len).*t_end;   
    one_path_times = min( v, A(lab_indices)');
    T_allow_times{2} = num2cell(one_path_times);
    
    %%% Above code does same as following for loop:    
    %%for ii = 1:len
    %%  v(ii) = min( v(ii) , A(lab{ii}(1), lab{ii}(2)));
    %%  end
    %%  % Now create fixed list of 1-path times and update T_allow_times{2}
    %%  T_allow_times{2} = num2cell(v);
    %%  one_path_times = v;

  
% Now for the other paths
for ii = 3:num_dims
  %  Set up the variables
  len = length(T_allow_times{ii});
  lab = T_labels{ii};
  v = ones(1,len).*t_end;
  % Get labels into matrix format
  lab_vector = cell2mat(lab);
  lab_array = reshape(lab_vector,ii,len)';
    %  Set up nested loop which computes allow times
    % column i of allow_times gives entry times of columns i, i+1 in lab_array
    allow_times = ones(len,ii-1).*t_end; 
    for jj = 1:ii-1
      lab_indices = sub2ind(size(A), lab_array(:,jj), lab_array(:,jj+1));
      allow_times(:,jj) = min(v', A(lab_indices));     
    end        
    %    Max of allow_times along columns is the entry time
    T_allow_times{ii}=num2cell(max(allow_times,[],2)');
  end
  
% Finally, sort by increasing allow times
  
  % Initialize
  T_sorted=cell(size(T_labels));
  % Populate and sort
  for kk=1:length(T_sorted)
    T_sorted{kk} = [T_labels{kk}', T_allow_times{kk}'];
    dummy_var = T_sorted{kk};
    [~,indx] = sort(cell2mat(dummy_var(:,2)));
    T_sorted{kk} = dummy_var(indx,:);
    end
    
% Output sorted labels and allow times

for kk = 1:length(T_sorted)
  T_labels{kk} = T_sorted{kk}(:,1)';
  T_allow_times{kk} = T_sorted{kk}(:,2)';  
  end


end