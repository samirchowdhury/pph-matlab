function [A,c] = pph_regularize(A)
%% This function takes an array and removes rows with repeated entries
% Inputs:
% A an array


% compare all-vs-all for each row using `bsxfun`
c = bsxfun( @eq, A, permute( A, [1 3 2] ) );
c = sum( c, 3 ); % count the number of matches each element has in the row
c = any( c > 1, 2 ); % indicates rows with repeated values - an element with more than one match
A = A( ~c, : );
end