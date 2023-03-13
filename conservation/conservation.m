%% energy conservation
% dx/vth is the typical time between collisions in a waterbag with vth
% the method to get the time of collision depends on the precision dt
% for dt<0.1 dx/vth the energy changes ~ 10^-7 in relative terms
% for dt~dx/vth enery starts to deviate significantly
% maximum relative error in energy seems to grow exponentially with dt for
% this specific setup

nparts = 1000;
dtfaclst = linspace(0.1,1.3,20);
err = 0*dtfaclst;
for i=1:numel(dtfaclst)
    disp(i)
    dtfac = dtfaclst(i);
    err(i) = conservation_aux(nparts,dtfac);
end

%% plot
semilogy(dtfaclst,err,'-b','LineWidth',2)
%% style
fnt = 24;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([2 1 1]) % 1.62
xlabel('$\kappa = dt \ vth/dx$','FontSize', fnt, 'Interpreter','latex')
ylabel('max$_t| \frac{\varepsilon(t)-\varepsilon(0)}{\varepsilon(0)}|$','FontSize', fnt, 'Interpreter','latex')
t=title('Maximum relative error of total energy','FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize';
xlim([min(dtfaclst),max(dtfaclst)])
yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
%% save plot
%print(gcf,'conservation.pdf','-dpdf','-r300')