%% get plot on waterbag distribution
addpath '..'
h = figure;
set(h,'PaperSize',[70 40]);

%% parameters
tmax = 45; % simulation time
dt = 0.1/50; % time step
tdim = round(tmax/dt); % number of iterations
nbins = 60; % number of bins for histogram
nparts = 40000; % number of particles
ndump = 1000; % time steps between dumps
ndumps = round(tdim/ndump); % number of dumps
boxL = 40; % plasma length
vfl = 0.4; % fluid velocity
vth = -0.1; % thermal velocity

% initial conditions
if vth>=0
    % maxwellian
    vi0 = vfl+normrnd(0,vth,1,nparts);
else
    % waterbag
    vi0 = vfl+vth*(2*rand(1,nparts)-1);
end

[histN,histedges] = histcounts(vi0,linspace(vfl-3*abs(vth),vfl+3*abs(vth),nbins));
[XS1,YS1] = histline(histedges, histN);
subplot(1,2,1)
plot(XS1,trapz(XS1,YS1)*1/(2*vth)*(heaviside(XS1-vfl+vth)-heaviside(XS1-vfl-vth)),'-b','LineWidth',2)
hold on
plot(XS1,YS1,'.k','MarkerSize',12)
%% style
fnt = 34;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([1 1 1]) % 1.62
xlabel('$v[c]$','FontSize', fnt, 'Interpreter','latex')
ylabel('$dN/dv [\textrm{arb. u.}]$','FontSize', fnt, 'Interpreter','latex')
t=title('Waterbag','FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize';
xlim([min(XS1),max(XS1)])
ylim([0,1.2*max(YS1)])

vth = -vth;
if vth>=0
    % maxwellian
    vi0 = vfl+normrnd(0,vth,1,nparts);
else
    % waterbag
    vi0 = vfl+vth*(2*rand(1,nparts)-1);
end
[histN,histedges] = histcounts(vi0,linspace(vfl-3*abs(vth),vfl+3*abs(vth),nbins));
[XS2,YS2] = histline(histedges, histN);
subplot(1,2,2)
plt2=plot(XS2,trapz(XS2,YS2)*1/sqrt(2*pi*vth^2)*exp(-0.5*((XS2-vfl)/vth).^2),'-b','LineWidth',2);
hold on
plt1=plot(XS2,YS2,'.k','MarkerSize',12);
%% style
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([1 1 1]) % 1.62
xlabel('$v[c]$','FontSize', fnt, 'Interpreter','latex')
%ylabel('$dN/dv [\textrm{arb. u.}]$','FontSize', fnt, 'Interpreter','latex')
yticks([])
t=title('Maxwellian','FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize';
legend([plt1,plt2],{'Code','Theory'},'FontSize',fnt, 'Interpreter','latex','Location','NorthWest')
xlim([min(XS1),max(XS1)])
ylim([0,1.2*max(YS1)])

sgtitle('Velocity distributions: vfl=0.4, vth=0.1','FontSize', fnt, 'Interpreter','latex')
%% save plot
%print(gcf,'waterbag.pdf','-dpdf','-r300')