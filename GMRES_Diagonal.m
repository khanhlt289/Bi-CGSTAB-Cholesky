function [x, converged, iter_ktl, res_norms] = GMRES_Diagonal(a_val, a_row_ptr, a_col_idx, b, restart, res_tol, max_iter, m_diagonal_val, use_HH)
% Generalized Minimum Residual Method with restarting 
% Correspond to Algorithm 6.9, 6.10, 6.11 in Yousef Saad's "Iterative Methods for Sparse Linear System (2nd Edition)"

n = size(a_row_ptr,1) -1 ; m = 10;
	
% 	if (nargin < 3) 
% 		restart  = min(n, 10); 
% 	end
% 	if (nargin < 4)
% 		res_tol  = 1e-9;
% 	end
% 	if (nargin < 5) 
% 		max_iter = min(floor(n / restart), 10); 
% 	end
% 	if (nargin < 6)
% 		M = eye(n);
% 	end
% 	if (nargin < 7) 
% 		use_HH = 1;  % Use Householder by default
% 	end 

	x = zeros(n, 1);
    r = b - dotproduct(a_val, a_row_ptr, a_col_idx, x);
% 	r = b - A * x_val;
    r = r ./ m_diagonal_val;

%     r = UpperTriangularSystem_M(M, r)
	residual = norm(r);
	stop_res = residual * res_tol;
	out_iter = 0;
	iter_ktl = 1; s = max_iter * restart + 1;
	res_norms = zeros(s, 1);
	res_norms(iter_ktl) = residual;
	
	converged = 0;
	while ((out_iter < max_iter) && (residual > stop_res))
		if (use_HH == 1) 
			[W, H, beta] = GMRES_Householder_Diagonal(a_val, a_row_ptr, a_col_idx, r, m, m_diagonal_val);
			[y, resvec]  = UpperHessenLeastSquare(H, -beta); 
			z = zeros(n, 1);
			for j = restart : -1 : 1
				% z = Pj * (z + e_j * y(j))
				z(j) = z(j) + y(j);
				z = z - 2 * W(:, j) * (W(:, j)' * z);
			end
		else
			[V, H, beta] = Arnoldi_MGS_Diagonal(a_val, a_row_ptr, a_col_idx, r, restart, m_diagonal_val);
			[y, resvec]  = UpperHessenLeastSquare(H, beta); 
			z = V(:, 1 : restart) * y(1 : restart);
		end
		
		x = x + z;
		r = b - dotproduct(a_val, a_row_ptr, a_col_idx, x);
		r = r ./ m_diagonal_val;
		
		for j = 1 : restart
			iter_ktl = iter_ktl + 1;
			res_norms(iter_ktl) = resvec(j);
			residual = min(residual, resvec(j));
			if (residual < stop_res)
				break;
			end
		end
	end
	if (residual <= stop_res) converged = 1; end
	
	res_norms = res_norms(1 : iter_ktl);
end