clear
close all

addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti")
addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1\distmesh\distmesh")
addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1")

h = 0.05;
[p, t, I, Bordo] = triangular_mesh('circle', h);

[p, t_Ext, e] = NVestendimesh(p, t);
[ne,~]=size(e);
[n, ~] = size(t);
[n_v, ~] = size(p);
n_i = size(I, 2);
n_b = size(Bordo, 2);
int_f = zeros(n, 1);
g_tilde = zeros(ne, 1);
A = spalloc(ne, ne, 8*n_v);
B = spalloc(n, ne, 8*n_v);

f = @(x) 1;
g = @(x) 1;

u = @(x) -1/4*(x(:, 1).^2+x(:, 2).^2)+5/4;
grad_u = @(x) -1/2*x;

for tr=1:n
    t1 = t(tr, 1);
    t2 = t(tr, 2);
    t3 = t(tr, 3);
    V1 = p(t1, :);
    V2 = p(t2, :);
    V3 = p(t3, :);

    b_t = (V1+V2+V3)./3;

    Area(tr) = polyarea([V1(1), V2(1), V3(1)], [V1(2), V2(2), V3(2)]);
    
    for iloc=1:3
        normale_i = t_Ext(tr, 6+iloc);
        lato_i = t_Ext(tr, 3+iloc);
        lunghezza_i = abs(p(e(lato_i, 1)) - p(e(lato_i, 2)));
        index_V_i = setdiff(t_Ext(tr, 1:3), e(lato_i, :));
        V_i = p(index_V_i, :);

        B(tr, lato_i) = B(tr, lato_i) - normale_i*lunghezza_i;

        for jloc=1:3
            jglob = t(tr, jloc);
            normale_j = t_Ext(tr, 6+jloc);
            lato_j = t_Ext(tr, 3+jloc);
            lunghezza_j = abs(p(e(lato_j, 1)) - p(e(lato_j, 2)));
            index_V_j = setdiff(t_Ext(tr, 1:3), e(lato_j, :));
            V_j = p(index_V_j, :);
           
            A(lato_i, lato_j) = A(lato_i, lato_j) + 1/(4*Area(tr))*(normale_i*normale_j...
                *lunghezza_i*lunghezza_j)*(b_t-V_i)*(b_t-V_j)';

 %           g_tilde(lato_j) = g_tilde(lato_j) + normale_i*normale_j*lunghezza_i...
 %           *lunghezza_j/2/Area(tr)*g(b_t)*(b_t-V_i)*(b_t-V_j)';
            g_tilde(lato_j) = g_tilde(lato_j) + g(b_t)*Area(tr);
       end
    end
    Bar(tr,:) = b_t;
    int_f(tr) = int_f(tr)-f(b_t)*Area(tr);    
end

AB = [A B'; B zeros(size(B,1), size(B,1))];
FG = [int_f; g_tilde];
UP = AB\FG;
U = UP(1:size(B,2),1);
P = UP(size(B,2)+1:end);

figure(2)
subplot(1,2,1)
plot3(Bar(:,1), Bar(:,2), P, '.')
xlabel('x')
ylabel('y')
zlabel('z')
title('Soluzione approssimata')
subplot(1,2,2)
plot3(Bar(:,1), Bar(:,2), u(Bar), '.')
xlabel('x')
ylabel('y')
zlabel('z')
title('Soluzione esatta')