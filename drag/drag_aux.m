function [tdump, vitfast, vitfastth, dx] = drag_aux(v,vth)
%% drag on fast sheet
% follow fast sheet index to get exact velocity (instead of Dawsons
% inaccurate method)

%% parameters
tmax = 30; % simulation time
dt = 0.1/50; % time step
tdim = round(tmax/dt); % number of iterations
nparts = 2000; % number of particles
ndump = 1000; % time steps between dumps
ndumps = round(tdim/ndump); % number of dumps
boxL = 40; % plasma length
vfl = 0.0; % fluid velocity

% initial conditions
if vth>=0
    % maxwellian
    vi0 = vfl+normrnd(0,vth,1,nparts);
else
    % waterbag
    vi0 = vfl+vth*(2*rand(1,nparts)-1);
end
% apply velocity of first sheet
fastindx = 1;
vi0(fastindx) = v;
vitfast = zeros(ndumps,1); % fast sheet velocity in time
vitfast(1) = vi0(fastindx);
xi0 = zeros(1,nparts); % initial sheet position
xieq = linspace(0,boxL,nparts); % sheet equilibrium position
dx = xieq(2)-xieq(1);
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
xlst = linspace(0,boxL,5*nparts); % x array for calculating field
Elst = zeros(ndumps,numel(xlst)); % field in xlst
Eilst = zeros(ndumps,numel(xlst)); % field @ sheet position

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
            % update fast sheet index?
            if (fastindx==i)
                fastindx = i+1;
            elseif(fastindx==i+1)
                fastindx = i;
            end
        end
    end
    % update previous position and velocity
    xipre = xi;
    
    % diagnostics
    if (mod(n,ndump)==0)
        tdump(idump) = t;
        vitfast(idump) = vi(fastindx);
        xit(idump,:) = xi;
        vit(idump,:) = vi;
        ene(idump) = sum(vi.^2 + (xi-xieq).^2);
        mom(idump) = sum(vi);
        [Elst(idump,:),Eilst(idump,:)] = getfields(xi,xlst,dx);
        idump = idump + 1;
    end
end

% theoretical drag on fast sheet
vitfastth = v-tdump*dx/2;