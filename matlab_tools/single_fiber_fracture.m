clc;clear;

D0 = 1; % diameter of fiber
a_by_c = 0.6;

if a_by_c == 0.6
    a_by_D = [0.05; 0.125; 0.2; 0.275; 0.35];
    F_values = [1.107; 1.176; 1.316; 1.565; 1.835]; % Table 3 raju and Newman
    F_fit = fit(a_by_D, F_values, 'poly2');
end

Q = 1 + 1.464*(a_by_c)^1.65;

a_range = linspace(0.05 , 0.35* D0);
K_range = linspace(0 , 1.6); % MPa.m^0.5

F_range = a_range / D0;
F_range = F_fit(F_range);

[a,K] = meshgrid(a_range, K_range);
[F1,dum2] = meshgrid(F_range, K_range);

a_val = realsqrt(a);
factor = sqrt(Q/pi)./F1;
stress = K ./ a_val;
stress = stress .* factor;


sol_a = linspace(0.05 , 0.35);
sol_K = F_range' .* sqrt(pi*sol_a/Q);


figure(1)
[C,h] = contourf(a,K,stress,250);
set(h,'LineColor','none')
xlabel('a/D_0');
ylabel('K/(S_0\surdD_0)');
colorbar;
colormap(jet);
caxis([0,1])
hold on;
plot(sol_a, sol_K, 'linewidth',1, 'color','black');
hold off;

figure(2)
plot(sol_a, sol_K, 'linewidth',2, 'color','blue');
xlabel('a/D_0');
ylabel('K/(S_0\surdD_0)');

