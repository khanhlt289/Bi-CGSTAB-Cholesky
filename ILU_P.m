% Chuong trinh phan tich ma tran A theo phuong phap ILU_P
% Giai phuong trinh dai so tuyen tinh Ax = b
% Nguoi lap trinh: Luu Truong Khanh
% Ngay thuc hien: 8-3-2021
% Ket qua: Thanh cong!
function [a] = ILU_P(a, p)
% ma tran U la phan tu cua ma tran a
% ma tran L la he so cua bien doi ma tran a
n = size(a,1); p = 13;
lev = zeros(n,n);
for i = 1:n
    for j =1:n
        if a(i,j) ~= 0
            lev(i,j) = 0;
        else
            lev(i,j) = inf;
        end
    end
end

for i = 2:n
%     i;
%     i-1;
    for k = 1:i-1
%         k;
%         ik = [i,k,a(i,k)];
        levik = lev(i,k);
%         p;
        if lev(i,k) <= p
%             ika = [i,k,a(i,k)];
%             kka = [k,k,a(k,k)];
            a(i,k) = a(i,k) / a(k,k);
            for j = k+ 1 : n
%                 j;
%                 ija = [i,j,a(i,j)];
%                 ika = [i,k,a(i,k)];
%                 kja = [k,j,a(k,j)];
                a(i,j) = a(i,j) - a(i,k) * a(k,j);
                levij = lev(i,j) ;
                levik = lev(i,k);
                levkj = lev(k,j);
                levijk = lev(i,k) + lev(k,j) + 1;
                lev(i,j) = min(levij, levijk);
            end
        end
    end
    for j =1:n
%         lev(i,j);
        if lev(i,j) > p
%             ij = [i,j];
            a(i,j) = 0;
        end
    end
end
