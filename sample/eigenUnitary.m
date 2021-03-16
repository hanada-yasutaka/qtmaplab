clear all
if ismac
    addpath("/Users/hanada/OneDrive/Packages/qtmaplab/");
elseif isunix
    %addpath("/home/hanada/OneDrive/Packages/qtmaplab/");
    addpath("/nfs/qtmaplab/");
    addpath("/nfs/AdvanpixMCT-4.8.3.14440/");
end

dim = 50;
domain = [-pi pi;-pi pi];
basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);
sU = SplitUnitary(dim, domain, basis);
k = 1;
tau = 1;

matT = sH.matT(@funcT);
matV = sH.matV(@funcV);

matH = matT + matV;
[hevecs, hevalsmat] = eig(matH);
[hevals, sindex] = sort(real(diag(hevalsmat)));

hevecs = hevecs(:, sindex);
hstates= eigs2states(sH, hevecs, hevals);

matU = sU.mat_expVTV(@funcT, @funcV, tau);

[uevecs, uevals] = eig(matU);
sindex = sortindex(uevecs, hevecs);
uevecs = uevecs(:, sindex);
uevals = diag(uevals(sindex, sindex) );
%uevals = uevals(sindex);

ustates = eigs2states(sU, uevecs, uevals);

sample = 100;
tmax = 300;
twopi = 2*pi;
%q0 = (rand(1, sample)-0.5)*twopi;
p0 = (rand(1, sample)-0.5)*twopi*2;
q0 = rand(1, sample)*twopi;
%p0 = rand(1, sample)*twopi;
%linspace(-pi, pi, sample);
%p0 = zeros(1, sample);
trajq = [];
trajp = [];
for i=1:tmax
    [q0,p0] = fVTV(q0,p0, k,tau);
    q0 = q0 - floor((q0 - pi)/twopi)*twopi - twopi;
    %p0 = p0 - floor((p0 - pi)/twopi)*twopi - twopi;
    %q0 = q0 - floor(q0/twopi)*twopi;
    %p0 = p0 - floor(p0/twopi)*twopi;    
    trajq = [trajq, q0];
    trajp = [trajp, p0];
end    
%plot(q0, p0, '.')

for i=1:dim
    s = ustates(i); 
    
    tiledlayout(2,2)
    nexttile;
    plot(s.q, log10(abs2( s.qrep() )), '-')
    title(sprintf("norm=%f", norm(s.qrep()) ));
    
    nexttile;
    theta = linspace(-pi, pi, 100);
    z = exp(1j*theta);
    plot(real(z), imag(z), '-k');
    hold on 
    scatter(real(uevals), imag(uevals) );
    scatter(real(s.eigenvalue), imag(s.eigenvalue),100, 'filled');
    hold off
    nexttile;
    
    [x,y,z] = s.hsmrep();
    contour(x, y, z, 10, 'LineWidth', 3);
    hold on;
    d = scatter(trajq, trajp, 1, '.');%, 'MarkerSize', 1, 'Marker', 'o')
    %d = scatter(trajq, trajp+2*pi, 1, '.');%, 'MarkerSize', 1, 'Marker', 'o')    
    axis(reshape(domain.', 1, []));
    %axis([0 2*pi 0 2*pi])
    hold off
    
    nexttile
    plot(log10(abs2( s.prep() )), s.p, '-')    
    title(sprintf("norm=%f", norm(s.prep()) ))
    %e = s.eigenvalue;
    %str = sprintf("%d-th eigs,E_n=%f", i, e); 
    %title(str);
    %fprintf("press enter to the next:\n")
    fprintf("%d\n", i);
    %savefig(sprintf("test/test_%d.png", i));
    waitforbuttonpress
end

return



function [qq,pp] = fTV(q,p,k,tau)
    pp = p - k*dfuncV(q)*tau;
    qq = q + dfuncT(pp)*tau;
end

function [qq,pp] = fVT(q,p,k,tau)
    qq = q + dfuncT(p)*tau;
    pp = p - k*dfuncV(qq)*tau;
end

function [q,p] = fVTV(q,p,k,tau)
    p = p - k*dfuncV(q)*tau/2;
    q = q + dfuncT(p)*tau;
    p = p - k*dfuncV(q)*tau/2;
end


function [q,p] = fTVT(q,p,k,tau)
    q = q + dfuncT(p)*tau/2;
    p = p - k*dfuncV(q)*tau;
    q = q + dfuncT(p)*tau/2;    
end


function y = funcV(x)
    y = cos(x);
end

function y = funcT(x)
    y = x .^2 / 2;
end


function y = dfuncV(x)
    y = -sin(x);
end

function y = dfuncT(x)
    y = x ;
end
