% Chuong trinh phan tu tren duong cheo chinh cua ma tran tien dieu kien. 
% Nguoi lap trinh: Luu Truong Khanh
% Ngay cap nhat va sua chua: 8-3-2021
function [m_diagonal_val, u_ptr] = diagonal_element_5(a_val, a_row_ptr, a_col_idx)
    n = size(a_row_ptr,1) - 1;
    m_diagonal_val = zeros(n,1);
    u_ptr = zeros(n,1);
    for i = 1:n
        j1 = a_row_ptr(i);
        j2 = a_row_ptr(i+1) - 1;
        if j2 == i - 1
            u_ptr(i) = i;
            m_diagonal_val(i) = 0;
        end
        for j = j1: j2
            k = a_col_idx(j);
            if (k == i)
                u_ptr(i) = j;
                m_diagonal_val(i) = a_val(j);
            end
        end
    end
end