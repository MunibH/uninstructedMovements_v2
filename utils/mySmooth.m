function [out,kernsd] = mySmooth(x, N, varargin)
% operates on first dimension only
% gaussian kernel with window size N, std dev sigma
% varargin:
%   'reflect' - reflect N elements from beginning of time series across
%               0 index
%   'zeropad' - pad N zeros across 0 index

% returns:
% - out: filtered data
% - kernsd: std dev of gaussian kernel

if N == 0 || N == 1 % no smoothing
    out = x;
    kernsd = 0;
    return
end

if isrow(x) % if 1-d time series, ensure column vector
    x = x';
end

if nargin > 2
    % handle a boundary condition
    bctype = varargin{1};
    if strcmpi(bctype,'reflect')
        try
        x_filt = cat(1,x(1:N,:),x);
        catch
            'a'
        end
        trim = N + 1;
    elseif strcmpi(bctype,'zeropad')
        try
        x_filt = cat(1,zeros(N,1),x);
        catch
            'a'
        end
        trim = N + 1;
    elseif strcmpi(bctype,'none')
        x_filt = x;
        trim = 1;
    else
        warning('invalid boundary condition type - regular smoothing');
        x_filt = x;
        trim = 1;
    end
else
    % no bc handling
    x_filt = x;
    trim = 1;
end



Ncol = size(x_filt, 2);
Nel = size(x_filt, 1);

kern = gausswin(N);
kernsd = std(1:N);


kern(1:floor(numel(kern)/2)) = 0; %causal
kern = kern./sum(kern);

out = zeros(Nel, Ncol);
for j = 1:Ncol
    out(:, j) = conv(x_filt(:, j), kern, 'same');
end
out = out(trim:end,:);

end