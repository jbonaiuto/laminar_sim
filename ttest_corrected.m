function [tval,p] = ttest_corrected(x, varargin)

% Parse inputs
defaults = struct('correction', 0);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

m = 0;

% Process remaining arguments
tail = 0;    % code for two-sided
dim = find(size(x) ~= 1, 1);
if isempty(dim), dim = 1; end

nans = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end
df = max(samplesize - 1,0);
xmean = nanmean(x,dim);
varpop = nanvar(x,[],dim);
if params.correction==0
    params.correction=.01*max(varpop);
end
corrsdpop=sqrt(varpop+params.correction);
%ser = sdpop ./ sqrt(samplesize);
ser = corrsdpop ./ sqrt(samplesize);
tval = (xmean - m) ./ ser;
% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = 2 * tcdf(-abs(tval), df);
elseif tail == 1 % right one-tailed test
    p = tcdf(-tval, df);
elseif tail == -1 % left one-tailed test
    p = tcdf(tval, df);
end
