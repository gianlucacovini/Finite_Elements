function int = Q0(u, p, t)
% Prende u funzione vettoriale e gli altri termini da triangular_mesh
    [n, ~] = size(t); % n=numero di triangoli

    for i=1:n
        % calcolo il baricentro
        vert = t(i, :);
        i1 = vert(1); % indice del primo vertice
        i2 = vert(2);
        i3 = vert(3);

        V1 = p(i1, :); % coordinate del primo vertice
        V2 = p(i2, :);
        V3 = p(i3, :);
        B = (V1+V2+V3)/3;

        % valuto la u nel baricentro
        u_b(i) = u(B);
        
        % calcolo dell'area
        A(i) = 1/2*det([V1(1) V1(2) 1; V2(1) V2(2) 1; V3(1) V3(2) 1]);
        %A(i) = polyarea([V1(1) V2(1) V3(1)], [V1(2) V2(2) V3(2)]);
    end

    % calcolo l'integrale come somma
    int = A*u_b';
end