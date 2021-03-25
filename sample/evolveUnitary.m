clear all
AdvanpixMCT='AdvanpixMCT-4.8.3.14460/';
if ismac
    addpath("/Users/hanada/OneDrive/Packages/qtmaplab/");
    addpath(sprintf("/Users/hanada/Applications/%s",AdvanpixMCT));
elseif isunix
    [~, name] = system('hostname');
    if strcmp(strtrim(name), 'bohigas')
        addpath("/home/hanada/OneDrive/Packages/qtmaplab/");
        addpath(sprintf("/home/hanada/Applications/%s",AdvanpixMCT));        
    else
        addpath("/nfs/qtmaplab/");
        addpath(sprintf('/nfs/%s', AdvanpixMCT) );
    end
end

dim = 1200;
mp.Digits(150);
%domain = [-pi pi;-pi pi];
%domain = [0 2*pi;0 2*pi];
%domain = [-2*pi 2*pi;-pi pi];
%%domain = [0 2*pi;-pi pi];
%domain = [-100 100;-100 100];
domain = mp('[-4*pi 4*pi;-15 15]');
%domain = [-pi pi;-pi pi];


basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);
[sU, state] = SplitUnitary(dim, domain, basis);

k = 1;
%tau = 0.1
tau = mp('-0.3');

%matT = sH.matT(@funcT);
%matV = sH.matV(@funcV);
%matH = matT + matV;
%[hevecs, hevalsmat] = eig(matH);
%[hevals, sindex] = sort(real(diag(hevalsmat)));
%hevecs = hevecs(:, sindex);
%hstates= eigs2states(sH, hevecs, hevals);


sample = 100;
tmax = 300;
twopi = 2*pi;
q0 = linspace(-10, 10, sample);
p0 = zeros(1, sample);
%q0 = (rand(1, sample)-0.5)*twopi;
%p0 = (rand(1, sample)-0.5)*twopi*2;
%q0 = rand(1, sample)*twopi;
%p0 = rand(1, sample)*twopi;
%linspace(-pi, pi, sample);
%p0 = zeros(1, sample);
trajq = [];
trajp = [];

for i=1:tmax
    [q0,p0] = fVT(q0,p0, double(k), double(tau));
    %q0 = q0 - floor((q0 - pi)/twopi)*twopi - twopi;
    %p0 = p0 - floor((p0 - pi)/twopi)*twopi - twopi;
    %q0 = q0 - floor(q0/twopi)*twopi;
    %p0 = p0 - floor(p0/twopi)*twopi;    
    trajq = [trajq, q0];
    trajp = [trajp, p0];
end    


s = state.coherent(0, 0);
%s2 = state.coherent(mp('+pi'),mp('0'));

%s1 = state.coherent(-pi,0);
%s2 = state.coherent(pi,0);
%s = (s1 - s2)/sqrt(2);
%s.norm()
%s = state.coherent(mp('0'), mp('0'));
%s = hstates(1);
fig = figure();
axs = [subplot(2,2,1) subplot(2,2,2) subplot(2,2,3) subplot(2,2,4) ] ;

for i=0:1000
        
    if mod(i,100)==0
        tiledlayout(2,2)
        nexttile;
        %plot(sU.q, log10(abs2(vec)), '-');
        plot(s.q, log10(abs2( s.qrep() )), '-', 'LineWidth', 3)
        %title(sprintf("norm=%f", norm(vec) ));
        
        nexttile;
        nexttile;
        
        %s = eigs2states(sU, vec);
        vrange = [-2*pi 2*pi; -2*pi 2*pi];
        [x,y,z] = s.hsmrep('vrange', vrange, 'gridnum', 100);
        levels = linspace(-28, 1, 20);
        contourf(x, y, log10(z), levels);
        hold on;
        d = scatter(trajq, trajp, 0.1, '.');%, 'MarkerSize', 1, 'Marker', 'o')
        %d = scatter(trajq+2*pi, trajp, 0.1, '.');%, 'MarkerSize', 1, 'Marker', 'o')
        axis(reshape(vrange.', 1, []));
        hold off
        
        nexttile
        %plot(log10(abs2( ifft(vec) )), sU.p, '-');
        plot(log10(abs2( s.prep() )), s.p, '-', 'LineWidth', 3)
        %title(sprintf("norm=%f", norm(ifft(vec)) ))
        %e = s.eigenvalue;
        %str = sprintf("%d-th eigs,E_n=%f", i, e);
        %title(str);
        %fprintf("press enter to the next:\n")
        fprintf("%d\n", i);
        %savefig(sprintf("test/test_%d.png", i));
        waitforbuttonpress
        
        %vec = op(vec);
        %disp(evolve_op)
        %vec = evolve_op(vec)
        %disp(vec);
    end
    s = sU.expVTevolve(s, @funcT, @funcV, tau);
    %s = vec2state(sU, s);
    
end




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
    y = x.^2 /2 - 2 * cos(x);
end

function y = funcT(x)
    y = x .^2 / 2;
end


function y = dfuncV(x)
    y = x + 2 * sin(x);
end

function y = dfuncT(x)
    y = x ;
end
