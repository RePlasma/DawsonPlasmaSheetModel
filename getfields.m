function [Elst,Eilst]= getfields(xi,xlst,dx)
% Electric fields: Elst @ x , Ei @ xi
Elst = zeros(1,numel(xlst));
Eilst = zeros(1,numel(xlst));
for i=1:numel(xlst)
    x = xlst(i);
    Elst(i) = x-dx*sum(xi <= x);
end
for i=1:numel(xi)
    x = xi(i);
    Eilst(i) = x-dx*sum(xi < x);
end