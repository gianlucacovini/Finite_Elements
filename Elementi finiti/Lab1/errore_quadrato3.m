clear
close all

addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1\distmesh\distmesh")

%% Integrazione di u=e^(x+y) sul quadrato (h corretto)

u = @(x) (exp(1).^(x(1)+x(2)));
I = (exp(1)-1).^2;

n = 10;
H = logspace(-1.5, -0.8, n);

for i=1:n
    h = H(i);
    [p,t,~,~]=triangular_mesh('square', h);    

    I0 = Q0(u, p, t);
    I1 = Q1(u, p, t);
    I2 = Q2(u, p, t);

    E0(i) = abs(I-I0);
    E1(i) = abs(I-I1);
    E2(i) = abs(I-I2);
    
    [n_tr, ~] = size(t);
    for tr=1:n_tr
        V1 = tr(1);
        V2 = tr(2);
        V3 = tr(3);
        L() = abs(V1-V2)
end

figure(1)
loglog(H_corr, E0, '-o', H_corr, E1, '-d', H_corr, E2, '-*', H_corr, H_corr.^2, H_corr, H_corr.^3/1e5)
legend("Errore Q0", "Errore Q1", "Errore Q2", "y=h^2")
title("Andamento dell'errore di e^{(x+y)} sul quadrato (h corretto)")