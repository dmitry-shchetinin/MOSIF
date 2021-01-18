function A = mosif_set_entries_to_ones(A)
%% set values of nonzero entries of a sparse matrix to one
    [m, n] = size(A);
    [i, j, ~] = find(A);
    A = sparse(i, j, ones(length(i), 1), m, n);
end