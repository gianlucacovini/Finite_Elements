clear
close all

addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1\distmesh\distmesh")

%% Integrazione di u=1 sul cerchio unitario

fd=@(p) sqrt(sum(p.^2,2))-1;

u = @(x, y) 1;
I = pi;

n = 20;
h = logspace(-1.5, -0.8, n);

for i=1:n
    [p,t]=distmesh2d(fd,@huniform,h(i),[-1,-1;1,1],[]);    

    I0 = Q0(u, p, t);
    I1 = Q1(u, p, t);
    I2 = Q2(u, p, t);

    E0(i) = abs(I-I0);
    E1(i) = abs(I-I1);
    E2(i) = abs(I-I2);
end

figure(1)
loglog(h, E0, '-o', h, E1, '-d', h, E2, '*', h, h.^2)
legend("Errore Q0", "Errore Q1", "Errore Q2", "y=h^2")
title("Andamento dell'errore di u=1 sul cerchio unitario")