% Chuong trinh thu nghiem giai he phuong trinh dai so tuyen tinh Ax = b
% bang phuong phap Bi_CGSTAB co tien dieu kien phai K 
% Tien dieu kien phai K la 1 phan tich Cholesky thu nghiem Chol_trial.
% Ngay thu nghiem: 9-2-2021
% Ngay sua chua: 22-2-2021 (lan 1)
% Nguoi lap trinh va thu nghiem: Luu Truong Khanh
% Muc dich: Hoc tap.
% Ket qua: Thu nghiem da thanh cong! Nghiem xap xi da giai duoc.
% Nguoi kiem tra: Luu Truong Khanh. Ngay kiem tra: 10-2-2021.
% Phat trien phan mem nang cao: Luu Truong Khanh
% Ket qua phat trien va sua chua (lan 1): Thanh cong!
% Chay thu nghiem voi chuong trinh san UBoot 12x12m. Ket qua: Thanh cong!
% Chay thu nghiem voi chuong trinh tuong 2 lop. Ket qua: ...
% =====================================================================
% Chon A, b
clear all; clc;
% Lap he phuong trinh dai so tuyen tinh Ax = b

[a_val, a_row_ptr, a_col_idx] = sparse2csr(A,1);
% Chon tien dieu kien phai K
l = Chol_trial(A);
K = l * l';
[k_val, k_row_ptr, k_col_idx] = sparse2csr(K,1);
% Khoi tao gia tri x_0 ban dau
x_0 = zeros(n,1);
% r_0 = b - A* x_0;
r_0 = b - dotproduct(a_val, a_row_ptr, a_col_idx, x_0);
r_bar = r_0;
i = 1;
tol = 1e-10;
rho_0 = r_bar' * r_0;
p_1 = r_0;
% Giai he phuong trinh K*p_hat = p_1 theo phuong phap GMRES (tam thoi dung)
M = eye(n,n); %res_tol = 1e-10; max_iter = 1000;  use_HH = 0; restart = 10;
res_tol = 1e-10; max_iter = 1000;  use_HH = 0; restart = 10;
m_diagonal = ones(n,1); res_tol = 1e-10; max_iter = 1000;  use_HH = 0; restart = 10;
% [p_hat, ~, ~, ~] = GMRES(K, p_1, restart, res_tol, max_iter, M, use_HH);
[p_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, p_1, restart, res_tol, max_iter, m_diagonal, use_HH);

% v_1 = A * p_hat;
v_1 = dotproduct(a_val, a_row_ptr, a_col_idx, p_hat);
alpha_1 = rho_0 / (r_bar' * v_1);
s = r_0 - alpha_1 * v_1;
% [s_hat, ~, ~, ~] = GMRES(K, s, restart, res_tol, max_iter, M, use_HH);
[s_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, s, restart, res_tol, max_iter, m_diagonal, use_HH);

% t = A * s_hat;
t = dotproduct(a_val, a_row_ptr, a_col_idx, s_hat);

omega_1 = t' *s / (t' * t);
x_1 =  x_0 + alpha_1 * p_hat + omega_1 * s_hat;
r_1 = s - omega_1 * t;
r_bar = r_1;
i = i + 1;
res = x_1 - x_0;
norm_res = norm(res,2);
while ((i < max_iter) && (norm_res > tol))
    rho_1 = r_bar' * r_0;
    beta_1 = (rho_1 / rho_0) * (alpha_1 / omega_1);
    p_2 = r_1 + beta_1 * (p_1 - omega_1 * v_1);
%     [p_hat, ~, ~, ~] = GMRES(K, p_2, restart, res_tol, max_iter, M, use_HH);
    [p_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, p_2, restart, res_tol, max_iter, m_diagonal, use_HH);
    
%     v_2 = A * p_hat;
    v_2 = dotproduct(a_val, a_row_ptr, a_col_idx, p_hat);
    
    alpha_2 = rho_1 / (r_bar' * v_2);
    s = r_1 - alpha_2 * v_2;
%     [s_hat, ~, ~, ~] = GMRES(K, s, restart, res_tol, max_iter, M, use_HH);
    [s_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, s, restart, res_tol, max_iter, m_diagonal, use_HH);

%     t = A * s_hat;
    t = dotproduct(a_val, a_row_ptr, a_col_idx, s_hat);
    
    omega_2 = t' * s / (t' * t);
    x_2 = x_1 + alpha_2 * p_hat + omega_2 * s_hat;
    r_2 = s - omega_2 * t;
    r_1 = r_2;
    x_1 = x_2;
    p_1 = p_2;
    rho_0 = rho_1;
    alpha_1 = alpha_2;
    omega_1 = omega_2;
    v_1 = v_2;
    r_bar = r_1;
    i = i + 1;
    res = x_2 - x_1;
    norm_res = norm(res,2);
end
rr = A * x_2 - b;
norm(rr,2)