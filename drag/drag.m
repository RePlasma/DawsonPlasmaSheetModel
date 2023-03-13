%% Dawson1962 §C.1. Drag on a Fast Sheet
% comparison with theory
addpath '..'

vth = 0;
vlst = linspace(0.1,0.9,5);
for i=1:numel(vlst)
    v = vlst(i);
    disp(i)
    [tdump, vitfast, vitfastth, dx] = drag_aux(v,vth);
    plt1=plot(tdump,vitfast,'ok');
    hold on
    plt2=plot(tdump,vitfastth,'-b');
    hold on
end

%% style
fnt = 24;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([2 1 1]) % 1.62
xlabel('$t[\omega_p^{-1}]$','FontSize', fnt, 'Interpreter','latex')
ylabel('$v_{\textrm{sheet}}[c]$','FontSize', fnt, 'Interpreter','latex')
t=title(['Drag on a Fast Sheet: $\delta=$' sprintf('%.3g ',dx) '$c/\omega_p$'],'FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize';
legend([plt1,plt2],{'Code','Theory'},'FontSize',fnt, 'Interpreter','latex','Location','NorthEast')
xlim([0,max(tdump)])
ylim([0,1.1])
%% save plot
%print(gcf,'drag.pdf','-dpdf','-r300')