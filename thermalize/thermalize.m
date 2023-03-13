%% start system with a waterbag distribution, check evolution towards a maxwellian
% warning: this script takes a significant time to run

nparts = 1000*2;
boxL = 10; % plasma length
xieq = linspace(0,boxL,nparts); % sheet equilibrium position
dx = xieq(2)-xieq(1);

%% velocities
vfl = 0.0; % fluid velocity
vth = -0.3; % thermal velocity
% initial conditions
if vth>=0
    % maxwellian
    vi0 = vfl+normrnd(0,vth,1,nparts);
else
    % waterbag
    vi0 = vfl+vth*(2*rand(1,nparts)-1);
end

%% parameters
tmax = 12000; % simulation time
nbins = 40;
dt = 0.4*dx/abs(vth); % time step
tdim = round(tmax/dt); % number of iterations
ndump = 100; % time steps between dumps
ndumps = round(tdim/ndump)+1; % number of dumps

xi0 = zeros(1,nparts); % initial sheet position
% update initial conditions
xi = xi0 + xieq;
vi = vi0;
% safeguard previous position and velocity
xipre = xi;

%% diagnostics
ti0 = zeros(nparts,1); % sheet "phase"
tdump = zeros(ndumps,1); % time of dumps
xit = zeros(ndumps,nparts); % sheet positions in time
vit = zeros(ndumps,nparts); % sheet velocities in time
ene = zeros(ndumps,1); % energy in time
mom = zeros(ndumps,1); % momentum in time
colls = zeros(tdim,1); % collisions in time

%% first diagnostics
xit(1,:) = xi;
vit(1,:) = vi;
ene(1) = sum(vi.^2 + (xi-xieq).^2);
mom(1) = sum(vi);

%% simulate
idump = 2; % dump iterator
for n=2:tdim
    t = (n-1)*dt;
    
    %% push particles to t
    for i=1:nparts
        xi(i) = xi0(i)*cos(t-ti0(i))+vi0(i)*sin(t-ti0(i)) + xieq(i);
        vi(i) = vi0(i)*cos(t-ti0(i))-xi0(i)*sin(t-ti0(i));
    end
  
    for i=1:nparts-1
        if( xi(i)>xi(i+1) )
            % compute dtc1
            dtc1 = dt*(xipre(i+1)-xipre(i))/(xipre(i+1)-xipre(i)+xi(i)-xi(i+1));
            % compute positions at t-dt+dtc1
            xi1 = xi0(i)*cos(t-ti0(i)+dtc1-dt)+vi0(i)*sin(t-ti0(i)+dtc1-dt)+xieq(i);
            xi2 = xi0(i+1)*cos(t-ti0(i+1)+dtc1-dt)+vi0(i+1)*sin(t-ti0(i+1)+dtc1-dt)+xieq(i+1);
            % compute dtc2
            dtc2 = (dt-dtc1)*(xi2-xi1)/(xi2-xi1+xi(i)-xi(i+1));
            % compute i<->i positions at t+dtc1+dtc2 "crossing"
            xi11 = xi0(i)*cos(t-ti0(i)+dtc1+dtc2-dt)+vi0(i)*sin(t-ti0(i)+dtc1+dtc2-dt)+xieq(i);
            xi22 = xi0(i+1)*cos(t-ti0(i+1)+dtc1+dtc2-dt)+vi0(i+1)*sin(t-ti0(i+1)+dtc1+dtc2-dt)+xieq(i+1);
            % compute i<->i velocities at t+dtc1+dtc2
            vi11 = vi0(i)*cos(t-ti0(i)+dtc1+dtc2-dt)-xi0(i)*sin(t-ti0(i)+dtc1+dtc2-dt);
            vi22 = vi0(i+1)*cos(t-ti0(i+1)+dtc1+dtc2-dt)-xi0(i+1)*sin(t-ti0(i+1)+dtc1+dtc2-dt);
            % update t0
            ti0(i) = t+dtc1+dtc2-dt;
            ti0(i+1) = t+dtc1+dtc2-dt;
            % update i<->j swaped "initial conditions"
            xi0(i) = xi22-xieq(i);
            vi0(i) = vi22;
            xi0(i+1) = xi11-xieq(i+1);
            vi0(i+1) = vi11;
            % correct position and velocity at t
            xi(i) = xi0(i)*cos(t-ti0(i))+vi0(i)*sin(t-ti0(i)) + xieq(i);
            vi(i) = vi0(i)*cos(t-ti0(i))-xi0(i)*sin(t-ti0(i));
            xi(i+1) = xi0(i+1)*cos(t-ti0(i+1))+vi0(i+1)*sin(t-ti0(i+1)) + xieq(i+1);
            vi(i+1) = vi0(i+1)*cos(t-ti0(i+1))-xi0(i+1)*sin(t-ti0(i+1));
            % count collisions
            colls(n) = colls(n)+1; % collisions
        end
    end
    % update previous position and velocity
    xipre = xi;
    
    % diagnostics
    if (mod(n,ndump)==0)
        disp(t)
        tdump(idump) = t;
        xit(idump,:) = xi;
        vit(idump,:) = vi;
        ene(idump) = sum(vi.^2 + (xi-xieq).^2);
        mom(idump) = sum(vi);
        idump = idump + 1;
    end
end
% last dump
tdump(idump) = t;
xit(idump,:) = xi;
vit(idump,:) = vi;
ene(idump) = sum(vi.^2 + (xi-xieq).^2);
mom(idump) = sum(vi);
idump = idump + 1;

%{
plot(xit(1,:),vit(1,:),'.b')
hold on
plot(xit(end,:),vit(end,:),'.r')
%}

vthth = abs(vth)/2;
%% plots
% t=0
[histN,histedges] = histcounts(abs(vit(1,:)),linspace(0,vfl+2*abs(vth),nbins));
[XS1,YS1] = histline(histedges, histN);
plt1=plot(XS1,2*trapz(XS1,YS1)*1/(2*vth)*(heaviside(XS1-vfl+vth)-heaviside(XS1-vfl-vth)),'-b','LineWidth',2);
hold on
plt2=plot(XS1,YS1,'ob','LineWidth',2);
hold on
% t=tmax
[histN,histedges] = histcounts(abs(vit(end,:)),linspace(0,vfl+2*abs(vth),nbins));
[XS2,YS2] = histline(histedges, histN);
plt3=plot(XS2,2*trapz(XS2,YS2)*1/sqrt(2*pi*vthth^2)*exp(-0.5*((XS2-vfl)/vthth).^2),'-r','LineWidth',2);
hold on
plt4=plot(XS2,YS2,'or','LineWidth',2);
%% style
fnt = 24;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([2 1 1]) % 1.62
xlabel('$v[c]$','FontSize', fnt, 'Interpreter','latex')
ylabel('$dN/dv [\textrm{arb. u.}]$','FontSize', fnt, 'Interpreter','latex')
t=title('Thermalization of a Waterbag','FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize';
legend([plt1,plt3],{'$t=0$','$t=t$max'},'FontSize',fnt, 'Interpreter','latex','Location','NorthEast')
xlim([min(XS2),max(XS2)])
%% save plot
%print(gcf,'thermalize.pdf','-dpdf','-r300')

%% energy conservation?
ene(1)
ene(end)