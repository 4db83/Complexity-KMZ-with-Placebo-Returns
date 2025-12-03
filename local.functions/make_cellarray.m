function matcell_out = make_cellarray(mat_in)
% make nR by 1 cell array
  [nR,nC] = size(mat_in);
  matcell_out = mat2cell(mat_in, ones(nR,1), nC);
end