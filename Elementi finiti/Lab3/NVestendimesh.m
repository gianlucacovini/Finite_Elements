function [p, t_Ext, e] = NVestendimesh(p, t)

% Usare comando unique riordinando prima i lati
[n, ~] = size(t);
e = [];
Ext_lati = [];
Ext_norm = [];

for tr=1:n
    t1 = t(tr, 1);
    t2 = t(tr, 2);
    t3 = t(tr, 3);

    L(1, :) = [t2 t3];
    L(2, :) = [t3 t1];
    L(3, :) = [t1 t2];

    for iloc=1:3
        if all(size(e)==[0 0])
            ind_lato(iloc) = size(e, 1)+1;
            val_norm(iloc) = 1;
            e = [e; L(iloc, :)];            
        elseif any(ismember(e, [L(iloc, 2) L(iloc, 1)], "rows"))
            ind_lato(iloc) = find(all((e==[L(iloc, 2) L(iloc, 1)])'));
            val_norm(iloc) = -1;
        else
            ind_lato(iloc) = size(e, 1)+1;
            val_norm(iloc) = 1;
            e = [e; L(iloc, :)];
        end
    end
    Ext_lati = [Ext_lati; [ind_lato(1) ind_lato(2) ind_lato(3)]];
    Ext_norm = [Ext_norm; [val_norm(1) val_norm(2) val_norm(3)]];
end

t_Ext = [t Ext_lati Ext_norm];
end