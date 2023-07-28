
function uniqueMatrix = removeDuplicateRows(matrix)
% Removes duplicate rows from a matrix
% Input: matrix - a matrix with potentially duplicate rows
% Output: uniqueMatrix - a matrix with all duplicate rows removed

% Find unique rows using the 'unique' function
[uniqueRows, ~, rowIdx] = unique(matrix, 'rows');

% Determine the number of occurrences of each unique row using 'accumarray'
rowCounts = accumarray(rowIdx, 1);

% Select rows that occur only once
uniqueRows = uniqueRows(rowCounts == 1, :);

% Return the unique matrix
uniqueMatrix = uniqueRows;
end