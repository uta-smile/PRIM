function [ data, fdata ] = load_db( db_id , pars)
%
%
%
    switch db_id
        case 1
            db = load('rawdata_brain');
            fdata = db.raw_data;
        case 2
            db = load('brain_8ch.mat');
            fdata = db.DATA;
        case 3
            db = load('phantom.mat');
            fdata = db.DATA;
        case 4
            db = load('cardiac_perf_R8.mat');
            fdata = squeeze(db.kdata(:, :, 1, :));
        case 5
            db = load('knee_pMRI_s1.mat');
            fdata = db.fdata;
        otherwise
            error('Unknown database_id');
    end
    
    if nargin > 1
        if isfield(pars, 'resize')
            [~, ~, T] = size(fdata);
            data = ifft2c(fdata);
            for t = 1:T
                rdata(:, :, t) = imresize(data(:, :, t), pars.resize);
            end
            fdata = fft2c(rdata);
        end
        if isfield(pars,'normalize')
            if pars.normalize
                m = max(abs(fdata(:)));
                fdata = fdata/m;
            end
        end
        if isfield(pars, 'noise')
            fdata = fdata + pars.noise * complex(randn(size(fdata)), randn(size(fdata)));
        end
    end
    data = ifft2c(fdata);

end

