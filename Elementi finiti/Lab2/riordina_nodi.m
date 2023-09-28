function E_n = riordina_nodi(t, B_I)
% I è il bordo di Neumann, B il bordo di Dirichlet
    n = size(t, 1);
    E_n = [];
    for tr=1:n
        node=t(tr, :);
        node_neumann = intersect(node, B_I);
        if size(node_neumann, 2) == 2
            E_n = [E_n; node_neumann];
        end
    end
end
