function colIndices = findColumnsBySum(matrix, sumValue)
% Finds the column indices of all columns in a matrix which add up to a specified value
% Input: matrix - the matrix to search
%        sumValue - the target sum value
% Output: colIndices - an array containing the column indices of columns that add up to sumValue

% Compute the sum of each column in the matrix
colSums = sum(matrix);

% Find the indices of columns that have a sum equal to sumValue
colIndices = find(colSums == sumValue);

end
