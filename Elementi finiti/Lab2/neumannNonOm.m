clear
close all

addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti")
addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1\distmesh\distmesh")
addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1")


f = @(x) exp(1)^(2*x(1)+x(2));
k = @(x) x(1).^2+x(2).^2+1;
%k = @(x) 1;
g = @(x) 1;
eta = @(x) 100;

h = 0.05;
[p, t, I, B] = triangular_mesh_modified('square', h);

tol=1e-3;
D = [];
N = [];
n_b = size(B, 2);
for i=1:n_b
    index = B(i);
    if p(index,1) <= (-1+tol)
        D = [D; index];
    else
        N = [N; index];
    end
end

[n, ~] = size(t);
[n_v, ~] = size(p);
int_f_tilde = zeros(n_v, 1);
int_neu_tilde = zeros(n_v, 1);
A_tilde = spalloc(n_v, n_v, 8*n_v);

E_n = riordina_nodi(t, N);

for tr=1:n
    t1 = t(tr, 1);
    t2 = t(tr, 2);
    t3 = t(tr, 3);
    V1 = p(t1, :);
    V2 = p(t2, :);
    V3 = p(t3, :);

    b_t = (V1+V2+V3)./3;

    Area(tr) = polyarea([V1(1), V2(1), V3(1)], [V1(2), V2(2), V3(2)]);

    L(1, :) = V3(:)-V2(:);
    L(2, :) = V1(:)-V3(:);
    L(3, :) = V2(:)-V1(:);
    for iloc=1:3
        iglob = t(tr, iloc); 

        for jloc=1:3
            jglob = t(tr, jloc);
            
            integral = dot(L(jloc, :), L(iloc, :))/(4*Area(tr)); 
            
            A_tilde(iglob, jglob) = A_tilde(iglob, jglob) + integral;
        end

        int_f_tilde(iglob) = int_f_tilde(iglob)+Area(tr)/3*f(b_t);
    end
end

for segm=E_n'
    int_neu_tilde(segm(1)) = int_neu_tilde(segm(1))+(segm(2)-segm(1))*eta((p(segm(1),:)+p(segm(2),:))/2);
end

n_i = size(I, 2);
n_d = size(D, 2);

for vert_ind=1:n_v
    if vert_ind <= n_i
        g_tilde(vert_ind) = 0;
    else
        g_tilde(vert_ind) = g(p(vert_ind, :));
    end
end

int_neu = int_neu_tilde(I);
A_II = A_tilde(I, I);
A_ID = A_tilde(I, D);
int_f = int_f_tilde(I);
int_f = int_f+int_neu;
g_D = g_tilde(D)';

u_h_int = A_II\(int_f-A_ID*g_D);
u_h = [u_h_int; g_D];

figure(2)
trimesh(t, p(:,1), p(:,2), u_h)
xlabel('x')
ylabel('y')
zlabel('z')
title('Soluzione approssimata')