addpath '..'
h = figure;
set(h,'PaperSize',[70 40]);

%%
vth = 0; v = 2;
nparts = 100;
boxLlst = v*linspace(nparts-0.4*nparts,nparts+0.4*nparts,11);
xlstdim = 50*nparts;
Emax = boxLlst*0;
vmax = boxLlst*0;
collsavg = boxLlst*0;
Elst = zeros(numel(boxLlst),xlstdim);
vilst = zeros(numel(boxLlst),nparts);
xilst = zeros(numel(boxLlst),nparts);

%%
for i=1:numel(boxLlst)
    disp(i)
    boxL = boxLlst(i); % plasma length
    [Elst(i,:),vilst(i,:),xilst(i,:),collsavg(i)] = breakcold_aux(v,vth,nparts,xlstdim,boxL);
    Emax(i) = max(Elst(i,:)-mean(Elst(i,:)));
    indxs = vilst(i,:)~=v;
    vmax(i) = max( vilst(i,indxs) );
end



%% plot
subplot(1,3,1)
dxlst = boxLlst/nparts;
plot(dxlst/v,vmax,'.k','MarkerSize',12)
hold on
plot(dxlst/v,dxlst,'-b','LineWidth',2)
hold on
plot([1,1],[0,max(vmax)],'--b','LineWidth',2)
%% style
fnt = 34;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([1 1 1]) % 1.62
xlabel('$dx/v_\phi[\omega_p^{-1}]$','FontSize', fnt, 'Interpreter','latex')
ylabel('max$(v)[c]$','FontSize', fnt, 'Interpreter','latex')
xlim([min(dxlst),1.01*max(dxlst)]/v)

subplot(1,3,2)
plot(dxlst/v,Emax,'.k','MarkerSize',12);
hold on
plot(dxlst/v,1.5*dxlst,'-b','LineWidth',2);
hold on
plot([1,1],[0,max(Emax)],'--b','LineWidth',2)
%% style
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([1 1 1]) % 1.62
xlabel('$dx/v_\phi[\omega_p^{-1}]$','FontSize', fnt, 'Interpreter','latex')
ylabel('max$(E)$','FontSize', fnt, 'Interpreter','latex')
xlim([min(dxlst),1.01*max(dxlst)]/v)

subplot(1,3,3)
plt1=plot(dxlst/v,collsavg,'.k','MarkerSize',12);
hold on
plt2=plot([1,1],[0,max(collsavg)],'--b','LineWidth',2);
%% style
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([1 1 1]) % 1.62
xlabel('$dx/v_\phi[\omega_p^{-1}]$','FontSize', fnt, 'Interpreter','latex')
ylabel('$\langle$ \#coll $\rangle$','FontSize', fnt, 'Interpreter','latex')
legend([plt1,plt2],{'Code','Theory'},'FontSize',0.7*fnt, 'Interpreter','latex','Location','NorthWest')
xlim([min(dxlst),1.01*max(dxlst)]/v)

sgtitle(['Transition to (cold) wavebreaking: \#part=' sprintf('%.3g ',nparts)],'FontSize', fnt, 'Interpreter','latex');

%% save plot
%print(gcf,'breakcold.pdf','-dpdf','-r300')