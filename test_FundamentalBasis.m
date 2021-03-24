addpath("/nfs/AdvanpixMCT-4.8.3.14440/");
%addpath("/nfs/ResearchData/2020-A3/matlab/qmapix/")

dps = 100;

mp.Digits(dps);
format longG;

dim = 100;
domain = [-pi pi;-2*pi 2*pi];if ismac
    addpath("/Users/hanada/OneDrive/Packages/qtmaplab/");
    addpath(sprintf("/home/hanada/Applications/%s",AdvanpixMCT));
elseif isunix
    [~, name] = system('hostname');
    if strcmp(strtrim(name), 'bohigas')
        addpath("/home/hanada/OneDrive/Packages/qtmaplab/");
        addpath(sprintf("/home/hanada/Applications/%s",AdvanpixMCT));        
    else
        addpath("/nfs/qtmaplab/");
        addpath("/nfs/AdvanpixMCT-4.8.3.14440/");
    end
end

%domain = mp('[-pi pi; -2*pi, 2*pi]');
basis = 'q';
H = SplitHamiltonian(dim, domain, basis);
U = SplitUnitary(dim, domain, basis);

matT = H.qBaseMatT(@(x) x .* x / 2 );
matV = H.matV(@(x) cos(x));


matH = matV + matT;
[evecs, evalsmat] = eig(matH);
[evals, sindex] = sort(real(diag(evalsmat)));
evecs = evecs(:,sindex);

vec = evecs(:, 1)
size(vec)
semilogy(H.q, abs2(vec))
return 
%states = U.eigen2states(evals, evecs);
%a = states(1).q2p();
%display(a)
%b = states(2);
%y = a + a;
%y.data
return 
%a = state(:,1);
%state.eigenvalue
%size(states)


%size(diag(evalsmat))
%utils.save_eigenvalues(H, diag(evalsmat), "test_dir");
%return 

%vec = evecs(:, 1)
%U.setop(@T, @V)
%set(U,'funcV', @(x) x * x)
%op = U.getop()
%vec = op(invec)
%U.evolveTV(invec)

%semilogy(H.q, abs2(vec))

%hold on
%for i = 1:10
%    vec = U.evolveTV(vec);
%    disp(i);
%    semilogy(H.q, abs2(vec),'-', 'Linewidth',2);
%end
%hold off
%H.test(T)
%return 
%norm(evecs(:,1))


%semilogy(H.p, abs2(evecs(:, 2)))

function y = T(x)
y = x .* x / 2;
end

function y = V(x)
y = cos(x);
end




