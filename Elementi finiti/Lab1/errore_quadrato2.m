clear
close all

addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1\distmesh\distmesh")

%% Integrazione di u=x^2+y^2 sul quadrato

u = @(x) (x(1)^2+x(2)^2);
I = 2/3;

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
end

figure(1)
loglog(H, E0, '-o', H, E1, '-d', H, E2, '-*', H, H.^2)
legend("Errore Q0", "Errore Q1", "Errore Q2", "y=h^2")
title("Andamento dell'errore di u=x^2+y^2 sul quadrato")