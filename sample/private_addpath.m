function private_addpath(AdvanpixMCT)
    if ismac
        addpath('/Users/hanada/OneDrive/Packages/qtmaplab/');
        addpath(sprintf('/Users/hanada/Applications/%s',AdvanpixMCT));
    elseif isunix
        [~, name] = system('hostname');
        if strcmp(strtrim(name), 'bohigas')
            addpath('/home/hanada/OneDrive/Packages/qtmaplab/');
            addpath(sprintf('/home/hanada/Applications/%s',AdvanpixMCT));
        else
            addpath('/nfs/qtmaplab/');
            addpath(sprintf('/nfs/%s', AdvanpixMCT));
        end
    end
end

