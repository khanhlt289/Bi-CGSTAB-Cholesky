function [V, H, beta] = Arnoldi_MGS_Diagonal(a_val, a_row_ptr, a_col_idx, r0, m, m_diagonal_val); %= Arnoldi_MGS(A, r0, m, M)
% Modified Gram-Schmidt for orthonormalizing Krylov subspace $K_{m}(A, r0)$
% V: the orthonormalized basis; H, beta: for GMRES

    beta = norm(r0, 2);
	n = size(a_row_ptr, 1)-1;
	V = zeros(n, m + 1);
	H = zeros(m + 1, m);
	V(:, 1) = r0 / beta;
    v_val = V(:,1);
	for j = 1 : m
% 		w = A * V(:, j);
        w = dotproduct(a_val,a_row_ptr,a_col_idx,v_val);
% 		w = M \ w;  % Left preconditioning
        
        % Giai he phuong trinh x = Mr
        w = w./ m_diagonal_val;
		for i = 1 : j
			H(i, j) = w' * V(:, i);
			w = w - H(i, j) * V(:, i);
		end
		H(j + 1, j) = norm(w, 2);
		if (H(j + 1, j) == 0) 
			break;
		end
		V(:, j + 1) = w / H(j + 1, j);
        v_val = V(:, j + 1);
	end
end