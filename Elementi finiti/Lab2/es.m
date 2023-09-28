clear
close all

addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti")
addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1\distmesh\distmesh")
addpath("C:\Users\gianl\OneDrive\Documenti\Elementi finiti\Lab1")


f = @(x) 1;

u = @(x) -1/4*(x(:, 1).^2+x(:, 2).^2)+1/4;
grad_u = @(x) -1/2*x;

num_index = 12;
H = logspace(-1.5, -0.8, num_index);
ErrL2_quad = zeros(num_index, 1);
SemiH1_quad = zeros(num_index, 1);
ErrL2_pm = zeros(num_index, 1);
SemiH1_pm = zeros(num_index, 1);

for index=num_index:-1:1
    h = H(index);
    [p, t, I, B] = triangular_mesh('circle', h);
    
    [n, ~] = size(t);
    [n_v, ~] = size(p);
    int_f_tilde = zeros(n_v, 1);
    A_tilde = spalloc(n_v, n_v, 8*n_v);
    
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
    
    A = A_tilde(I, I);
    int_f = int_f_tilde(I);
    
    n_i = size(I, 2);
    n_b = size(B, 2);

    u_h_int = A\int_f;
    u_h = [u_h_int; zeros(n_b,1)];
    u_sol = u(p);

    ErrInf(index) = norm(u_sol-u_h, 'inf');

    for tr=1:n
        t1 = t(tr, 1);
        t2 = t(tr, 2);
        t3 = t(tr, 3);
        V1 = p(t1, :);
        V2 = p(t2, :);
        V3 = p(t3, :);

        b_t = (V1+V2+V3)./3;
        vec_bt(tr,:) = b_t;

        m1 = (V2+V3)./2;
        m2 = (V3+V1)./2;
        m3 = (V1+V2)./2;

        L(1, :) = V3(:)-V2(:);
        L(2, :) = V1(:)-V3(:);
        L(3, :) = V2(:)-V1(:);

        U_bar(tr) = u(b_t);
        Uh_bar(tr) = (u_h(t1)+u_h(t2)+u_h(t3))/3;
        ErrL2_quad(index) = ErrL2_quad(index) + Area(tr)*(U_bar(tr)-Uh_bar(tr))^2;
        ErrL2(index) = ErrL2_quad(index);
        
        uh_med(1) = (u_h(t2)+u_h(t3))/2;
        uh_med(2) = (u_h(t3)+u_h(t1))/2;
        uh_med(3) = (u_h(t1)+u_h(t2))/2;
        ErrL2_pm(index) = ErrL2_pm(index) +...
            +Area(tr)*((u(m1)-uh_med(1))^2+...
            +(u(m2)-uh_med(2))^2+...
            +(u(m3)-uh_med(3))^2)/3;

        grad_u_bar = grad_u(b_t)';
        
        grad_uh_bar = [0; 0];
        for i=1:3
            ind = t(tr, i);
            grad_phi(ind, :) = ([0 -1; 1 0]*L(i, :)'/2/Area(tr))';
            grad_uh_bar = grad_uh_bar + u_h(ind)*grad_phi(ind, :)';
            grad_uh_pm = grad_uh_bar;
        end
        SemiH1_quad(index) = SemiH1_quad(index) +...
            +Area(tr)*sum((grad_u_bar-grad_uh_bar).^2);
        SemiH1(index) = sqrt(SemiH1_quad(index));
        ErrH1 = sqrt(ErrL2_quad + SemiH1_quad);

        grad_u_pm = [grad_u(m1); grad_u(m2); grad_u(m3)];

        SemiH1_pm(index) = SemiH1_pm(index) +...
            +Area(tr)*((grad_u_pm(1,1)-grad_uh_pm(1))^2+...
            +(grad_u_pm(2,1)-grad_uh_pm(1))^2+...
            +(grad_u_pm(3,1)-grad_uh_pm(1))^2)/3+...
            +Area(tr)*((grad_u_pm(1,2)-grad_uh_pm(2))^2+...
            +(grad_u_pm(2,2)-grad_uh_pm(2))^2+...
            +(grad_u_pm(3,2)-grad_uh_pm(2))^2)/3;
        ErrH1_pm(index) = sqrt(SemiH1_pm(index) + ErrL2_pm(index));
    end
    ErrL2 = sqrt(ErrL2_quad);
    ErrL2_pm = sqrt(ErrL2_pm);
end

figure(2)
subplot(2,1,1)
loglog(H, ErrInf, '-<', H, 1e-1*H.^2, '--', H, ErrL2, '-<', H, ErrH1, '-<', 'LineWidth', 1.5)
legend('Errore norma infinito', 'h^2', 'Errore norma L2', 'Errore norma H1', 'Location', 'best')
xlabel('h')
ylabel('Errore')
title('Errore in funzione di h (baricentri)')
subplot(2,1,2)
loglog(H, 1e-1*H.^2, '--', H, H, '--', H, ErrH1, '-<', H, ErrH1_pm, '-<', 'LineWidth', 1.5)
legend('h^2', 'h', 'Errore baricentri', 'Errore punti medi', 'Location', 'best')
xlabel('h')
ylabel('Errore')
title('Errore in funzione di h norma H1')

figure(3)
subplot(1,2,1)
trimesh(t, p(:,1), p(:,2), u_h)
xlabel('x')
ylabel('y')
zlabel('z')
title('Soluzione approssimata')
subplot(1,2,2)
trimesh(t, p(:,1), p(:,2), u(p))
xlabel('x')
ylabel('y')
zlabel('z')
title('Soluzione esatta')

figure(4)
trimesh(t, p(:,1), p(:,2), u(p)-u_h)
xlabel('x')
ylabel('y')
zlabel('z')
title('Errore nei vertici')