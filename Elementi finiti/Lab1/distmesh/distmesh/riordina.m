function [ p_riodinati , t_riordinati ] = riordina( df,p,t )
%riordina: orienta i triangoli e riordina i nodi 
% mettendo prima quelli interni (df <0)
% e in ultimo  quelli di bordo(df =0).

n_nodi = size(p,1);
n_triangoli = size (t,1);


% TO DO: orientamento dei triangoli in senso antiorario

% Riordinamento dei nodi
old2new = zeros(n_nodi,1); % nuova posizione dei nodi con indicizzazione vecchia
new2old = zeros(n_nodi,1); % vecchia posizione dei nodi con indicizzazione nuova
tol = 10^-3;
i_interno = 1;
i_bordo = n_nodi;
for i = 1:n_nodi
    if abs(df(p(i,:))) < tol
        old2new(i)=i_bordo;
        new2old(i_bordo)=i;
        i_bordo = i_bordo -1
    else
        old2new(i)=i_interno;
        new2old(i_interno)=i;
        i_interno = i_interno +1;
    end
end
p_riodinati = p(new2old,:);
t_riordinati = old2new(t);

