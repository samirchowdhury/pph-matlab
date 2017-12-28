function barcode = compute_pers(A,T_labels,T_allow_times)
% Main function for computing barcodes
% Inputs:  
% A is a network
% T_labels is a cell array. Entry i is a cell containing labels for (i-1) paths
% T_allow_times is a cell array of the same dimensions as T_labels. Each entry 
% of each cell is the allow time of the corresponding path.
  
% Initialize
barcode = cell(size(T_labels));

% Compute constants
max_path_dim = length(T_labels);

%% Initialize cells for storing marked paths, entry times, cycles, and cycle times
T_marked = T_allow_times;
T_entry_times = T_allow_times;
T_cycles = cell(size(T_marked));
T_cycle_times = cell(size(T_marked));
for ii = 1:length(T_marked)
  T_marked{ii} = zeros(size(T_marked{ii}));
  T_cycles{ii}=cell(size(T_marked{ii}));
  T_cycle_times{ii}=cell(size(T_marked{ii}));
  end


%% For 0-paths, boundary vector is 0. Mark and move on.
T_marked{1} = ones(size(T_marked{1}));

%% Proceed to p-paths for p \geq 1

for ii=2:max_path_dim %start ii from 2
  
  for jj=1:length(T_labels{ii})
    % Get current path and its allow time
    path = T_labels{ii}{jj};
    [boundary_vector, summands] = compute_boundary(path,T_labels{ii-1},T_marked{ii-1},T_allow_times{ii-1});
    [reduced_boundary_vector, max_summand_indx, max_summand_at] = basis_change(boundary_vector, T_cycles{ii-1}, T_labels{ii-1}, T_allow_times{ii-1});
    % Take max of path's allow time and summand allow time
    et = max(T_allow_times{ii}{jj}, max_summand_at);
    
    %% not in submitted algorithm, added this. After column reduction,
    % update entry time.
    T_entry_times{ii}{jj} = et;
    
    % if reduced_boundary_vector = 0, mark T_marked{ii}{jj}
    if (~any(reduced_boundary_vector))
      T_marked{ii}(jj) = 1;
    else 
      T_cycles{ii-1}{max_summand_indx} = reduced_boundary_vector;
      T_cycle_times{ii-1}{max_summand_indx} = et;
      %% populate pers (FILL IN)
      barcode{ii-1} = [barcode{ii-1}; T_entry_times{ii-1}{max_summand_indx}, et];
      
    end
    
  
  end
  end
  
% add infinite bars

for ii = 1:max_path_dim-1
  for jj=1:length(T_labels{ii})
    if (T_marked{ii}(jj) && isempty(T_cycle_times{ii}{jj}))
      barcode{ii} = [barcode{ii}; T_entry_times{ii}{jj}, inf];
    end
    end
    end

% discard trivial bars (can comment this out)
for ii = 1:max_path_dim - 1
  barcode{ii} = pph_regularize(barcode{ii});
  end

  
  end
  
function [boundary_vector, summands]  = compute_boundary(path,list_bd_paths,list_bd_marks, list_bd_times)
%  Inputs
% path is a p-path, e.g. [1 2 3]
% list_bd_paths is just T_labels{p-1}
% list_bd_marks is just T_marked{p-1}
% list_bd_times is just T_allow_times{p-1}
% This function computes a boundary vector and its summands, along with 
% the summand allow times and slot positions. This can then go into 
% the column reduction procedure. 

% Compute constants
len = length(path);
boundary_vector = zeros(size(list_bd_marks));
length_bd_path = length(list_bd_paths{1});
num_bd_paths = length(list_bd_paths);

% reshape list_bd_paths to get into array format
list_bd_paths = reshape(cell2mat(list_bd_paths),length_bd_path,num_bd_paths)';
%list_bd_paths =reshape(cell2mat(list_bd_paths),length_bd_path,len)';
% list_bd_paths =reshape(cell2mat(list_bd_paths),2,len)';


% Initialize summands cell array.
% Column 1 contains the summands.
% Column 2 contains the slot for the summand in T_labels{p-1}.
% Column 3 contains the allow time of the corresponding summand.
summands = cell(len,3);

%% okay now that we have separate function for summand info, some of this can be killed (DO THIS)

% calculate summands and build boundary_vector
for ii=1:len
  v = path;
  v(ii) = [];
  summands{ii,1} = v;
  [~,indx] = ismember(v,list_bd_paths,'rows');
  % if v is of the form [1 1], indx will be 0. if so, pass to next iteration
  if (~indx)
    continue
    end
  summands{ii,2}=indx;
  summands{ii,3}=list_bd_times{indx};
  % if unmarked, remove from boundary vector; i.e. skip to next iteration
  if (list_bd_marks(indx))
  % if marked, carry on
  boundary_vector(indx) = (-1)^(ii+1);
  end
  end
  end
  
  
  
function [reduced_boundary_vector, max_summand_indx, max_summand_at] = basis_change(boundary_vector, list_cycles, list_paths, list_times)
% Inputs
% boundary_vector is obtained from compute_boundary
% list_cycles is just T_cycles{p-1}
% list_times and list_paths are T_allow_times{p-1} and T_labels{p-1}

% Outputs
% max_summand_at is the allow time of the maximal summand (maximal in the sense of allow time)
% max_summand_indx is the slot number of the maximal summand
  
  reduced_boundary_vector = boundary_vector;
  [summands, max_summand_indx, max_summand_at] = compute_summand_info_from_vector(reduced_boundary_vector,list_paths,list_times);
  
  % perform column reduction
  while (any(reduced_boundary_vector))
    [~, max_summand_indx, max_summand_at] = compute_summand_info_from_vector(reduced_boundary_vector,list_paths,list_times);
    % check if this slot in T_cycles is empty
%    % summands{max_summand_at_indx,2} gives slot index in list_cycles
%    indx = summands{max_summand_at_indx,2};
    if (isempty(list_cycles{max_summand_indx}))
      break
    end
    % this point occurs after the break, so list_cycles slot is known 
    % to be nonempty. so kill this entry in summands (FILL IN)
    
    % list_cycles{max_summand_indx} should be a vector of the same length as reduced_boundary_vector,
    % and should have a nonzero entry (+-1) in position max_summand_indx
    reduction_vector = list_cycles{max_summand_indx};
    reduction_scalar = reduced_boundary_vector(max_summand_indx).*reduction_vector(max_summand_indx);
    reduced_boundary_vector = reduced_boundary_vector - (reduction_scalar.*reduction_vector);
    
    end
    
    
  end
  
  
function [summands, max_summand_indx, max_summand_at] = compute_summand_info_from_vector(path_vector,list_paths,list_times)
% Inputs
% path_vector is a p-path in vector form
% list_paths is a cell array of all the p-paths 
% list_times is a cell array containing allow times of all p-paths

% Given these inputs, compute a cell array with 3 columns.
% Column 1 contains the summands in path form 
% Column 2 contains the slot number for each summand in T_labels
% Column 3 contains the allow time of each summand 
% Also compute the allow times of each summand, and pick out the (unique) 
% greatest index of all the summands with max allow time

% Compute number of nonzero entries in path_vector
len =  sum(path_vector~=0);
% initialize
summands = cell(len,3);
mm=1;
for kk=1:length(path_vector)
  if (path_vector(kk)==0)
    continue
    end
  summands{mm,1} = list_paths{kk};
  summands{mm,2} = kk;
  summands{mm,3} = list_times{kk};
  mm=mm+1;
  end
  
  % extract summand allow times
  summand_times = cell2mat(summands(:,3));
  % find max allow time
  max_summand_at = max(summand_times);
  % find all summands with max allow time
  max_summand_at_indices = (summand_times == max_summand_at);
  % create a vector containing the slot numbers of summands attaining this allow time
  summands_vector = cell2mat(summands(:,2));
  % among these summands, pick the one with the largest slot number (rightmost entry in T_labels)
  max_summand_indx = max(summands_vector.*max_summand_at_indices);
  
  
  end
%  
%% extract summand allow times and summand with max allow time
%summand_times = cell2mat(summands(:,3))
%[~,max_summand_at]=max(summand_times);
%
%  
%% if boundary_vector is nonzero, perform column reduction 
%
%while (any(boundary_vector))
%  
%  end

