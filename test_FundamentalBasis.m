import FundamentalBasis
addpath("/nfs/AdvanpixMCT-4.8.0.14100/");

dps = 100;

mp.Digits(dps);
format longG;

dim = 100;
%domain = [-pi pi;-2*pi 2*pi];
domain = mp('[-pi pi; -2*pi, 2*pi]');
rep = 'q';
H = FundamentalBasis(dim, domain, rep)

matT = H.matT(@(x) x .* x / 2 );
matV = H.matV(@(x) cos(x));
disp(H.q)
matH = matV + matT;
[evecs, evalsmat] = eig(matH);
[evals, sindex] = sort(real(diag(evalsmat)));
evecs = evecs(:,sindex);
class(evecs)

semilogy(H.p, abs2(evecs(:, 2)))

function y = abs2(x)
y = abs(conj(x) .* x);
end



