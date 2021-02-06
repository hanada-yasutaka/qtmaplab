%addpath("/nfs/AdvanpixMCT-4.8.3.14440/");
addpath("/Users/hanada/Documents/MATLAB/qtmaplab/");

dim = 50;
domain = [-pi pi; -pi pi];
basis = 'q';
sH = SplitHamiltonian(dim, domain, basis);
matT = sH.matT(@(x) x .^2 /2);
matV = sH.matV(@(x) cos(x));
 
matH = matT + matV;
[evecs, evalsmat] = eig(matH);
[evals, sindex] = sort(real(diag(evalsmat)));

evecs = evecs(:,sindex);
states = sH.eigs2states(evals, evecs);

for i=1:dim
    s = states(i); 
    [x,y,z] = s.hsmrep();
    contour(x, y, z, 10);
    e = s.eigenvalue;
    str = sprintf("%d-th eigs,E_n=%f", i, e); 
    title(str);
    fprintf("press enter to the nesx:\n")
    waitforbuttonpress
end