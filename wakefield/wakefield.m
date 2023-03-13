addpath '..'

%%
vth = 0;
nparts = 2000;
xlstdim = 20*nparts;
vd = 40/nparts;
vlst = linspace(0.1,1,1);
Elst = zeros(numel(vlst),xlstdim);
vilst = zeros(numel(vlst),nparts);
xilst = zeros(numel(vlst),nparts);

%%
for i=1:numel(vlst)
    disp(i)
    v = vlst(i);
    [Edump1,Edump2,Edump3,Edump4,xlst,vilst(i,:),xilst(i,:)] = wakefield_aux(v,vth,nparts,xlstdim);
end

%% plot
amp = 2.1*max(Edump4-mean(Edump4));
plot(xlst,Edump1-mean(Edump1)+3*amp,'-b')
hold on
plot(xlst,Edump2-mean(Edump2)+2*amp,'-b')
hold on
plot(xlst,Edump3-mean(Edump3)+amp,'-b')
hold on
plot(xlst,Edump4-mean(Edump4),'-b')
%% style
fnt = 24;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([2 1 1]) % 1.62
title('Fast sheet')
xlabel('$x[c/\omega_p]$','FontSize', fnt, 'Interpreter','latex')
ylabel('$E$','FontSize', fnt, 'Interpreter','latex')
t=title('Wake of constant velocity sheet driver','FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize';
yticks([])
xlim([0,max(xlst)])
%% save plot
%print(gcf,'wakefield.pdf','-dpdf','-r300')