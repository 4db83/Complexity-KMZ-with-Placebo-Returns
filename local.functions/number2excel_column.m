function colLetter = number2excel_column(n)
%NUM2EXCELCOL Convert a positive integer to an Excel column name.
%   col = NUM2EXCELCOL(n) returns the Excel column label for positive
%   integer n (e.g. 1 -> 'A', 26 -> 'Z', 27 -> 'AA', 53 -> 'BA').

    if ~isscalar(n) || n <= 0 || n ~= floor(n)
        error('Input must be a single positive integer.');
    end

    letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M', ...
               'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];

    colLetter = '';
    while n > 0
        r = mod(n - 1, 26);
        colLetter = [letters(r + 1), colLetter]; %#ok<AGROW>
        n = floor((n - 1) / 26);
    end

% return string
colLetter = string(colLetter);

end
