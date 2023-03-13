addpath '..'

%%
vth = 0;
nparts = 2000;
xlstdim = 20*nparts;
vd = 40/nparts;
vlst = linspace(0.1,1,20);
Elst = zeros(numel(vlst),xlstdim);
vilst = zeros(numel(vlst),nparts);
xilst = zeros(numel(vlst),nparts);

%%
for i=1:numel(vlst)
    disp(i)
    v = vlst(i);
    [Elst(i,:),vilst(i,:),xilst(i,:)] = wavelength_aux(v,vth,nparts,xlstdim);
end

%% get wavelength for each sheet velocity
lbdlst = zeros(numel(vlst),1);
for i=1:numel(vlst)
    % remove fast sheet from phase space
    indxs = vilst(i,:)~=vlst(i);
    x = xilst(i,indxs);
    y = vilst(i,indxs);
    
    % https://www.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function
    yu = max(y);
    yl = min(y);
    yr = (yu-yl);                               % Range of ‘y’
    yz = y-yu+(yr/2);
    zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
    per = 2*mean(diff(zx));                     % Estimate period
    ym = mean(y);                               % Estimate offset
    fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
    fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
    s = fminsearch(fcn, [yr;  per;  -1;  ym]);                       % Minimise Least-Squares
    xp = linspace(min(x),max(x));
    
    % save wavelength
    lbdlst(i) = s(2);
    %plot(x,y,'b',  xp,fit(s,xp), 'r')
end

%% plot results
plt1=plot(vlst,lbdlst,'.k');
hold on
plt2=plot(vlst,2*pi*vlst,'-b');
%% style
fnt = 24;
ax = gca;
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.FontSize = fnt;
ax.TickLabelInterpreter = 'latex';
pbaspect([2 1 1]) % 1.62
title('Fast sheet')
xlabel('$v[c]$','FontSize', fnt, 'Interpreter','latex')
ylabel('$\lambda_p[c/\omega_p]$','FontSize', fnt, 'Interpreter','latex')
t=title('Wake of constant velocity sheet driver','FontSize', fnt, 'Interpreter','latex');
t.Units = 'Normalize'; 
%t.Position(1) = 0;
%t.HorizontalAlignment = 'left';
legend([plt1,plt2],{'Code','Theory'},'FontSize',fnt, 'Interpreter','latex','Location','NorthWest')
xlim([0,1])
ylim([0,2*pi])
%% save plot
%print(gcf,'wavelength.pdf','-dpdf','-r300')

%{
%% plot phase space
for i=1:numel(vlst)
    plot(xilst(i,:),vilst(i,:),'-')
    hold on
end
%}

%{
%% plot electric field
for i=1:numel(vlst)
    plot(Elst(i,:),'-')
    hold on
end
%}