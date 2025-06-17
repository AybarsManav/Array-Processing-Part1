function is_all_diagonal = check_diagonal_blocks(D_joint, tolerance)
%CHECK_DIAGONAL_BLOCKS Checks if each 2x2 block in a 2x(2m) matrix is diagonal.
%
%   is_all_diagonal = CHECK_DIAGONAL_BLOCKS(D_joint, tolerance)
%
%   Inputs:
%       D_joint   - A 2 x (2m) matrix composed of m horizontal 2x2 blocks.
%       tolerance - A scalar value for allowable off-diagonal error.
%
%   Output:
%       is_all_diagonal - Returns true if all blocks are diagonal within the tolerance.

    if nargin < 2
        tolerance = 1e-3;  % default value if not provided
    end

    [m_rows, m_cols] = size(D_joint);

    if m_rows ~= 2 || mod(m_cols, 2) ~= 0
        error('D_joint must be a 2 x (2m) matrix with 2x2 blocks.');
    end

    num_blocks = m_cols / 2;
    is_all_diagonal = true;

    for i = 1:num_blocks
        idx = (2*i - 1):(2*i);        % Column indices for i-th 2x2 block
        block = D_joint(:, idx);

        % Extract off-diagonal elements
        off_diagonals = [block(1,2), block(2,1)];

        if any(abs(off_diagonals) > tolerance)
            is_all_diagonal = false;
            fprintf('Block %d is NOT diagonal (tolerance %g):\n', i, tolerance);
            disp(block);
        end
    end

    if is_all_diagonal
        fprintf('All 2x2 blocks are diagonal within tolerance %.1e.\n', tolerance);
    else
        fprintf('Some 2x2 blocks are not diagonal within tolerance %.1e.\n', tolerance);
    end
end

