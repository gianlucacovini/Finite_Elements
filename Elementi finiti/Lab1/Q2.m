function int = Q2(u, p, t, I, B)
% Prende u funzione vettoriale e gli altri termini da triangular_mesh
    [n, ~] = size(t); % n=numero di triangoli

    % calcolo i punti medi dei lati e area per ogni triangolo
    for i=1:n        
        vert = t(i, :);
        i1 = vert(1); % indice del primo vertice
        i2 = vert(2);
        i3 = vert(3);

        V1 = p(i1, :); % coordinate del primo vertice
        V2 = p(i2, :);
        V3 = p(i3, :);

        m1 = (V1+V2)/2;
        m2 = (V2+V3)/2;
        m3 = (V3+V1)/2;
       
        % calcolo la valutazione di u
        u_val(i) = (u(m1)+u(m2)+u(m3))/3;

        A(i) = 1/2*det([V1(1) V1(2) 1; V2(1) V2(2) 1; V3(1) V3(2) 1]);
    end

    % calcolo l'integrale come somma
    int = A*u_val';
end