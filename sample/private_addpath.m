function private_addpath(AdvanpixMCT)
    if ismac
        addpath('/Users/hanada/Dropbox/Packages/qtmaplab/');
        addpath(sprintf('/Users/hanada/Applications/%s',AdvanpixMCT));
    elseif isunix
        [~, name] = system('hostname');
        else
            addpath('/nfs/Dropbox/Packages/qtmaplab/');
            addpath(sprintf('/nfs/%s', AdvanpixMCT));
        end
    end
end

