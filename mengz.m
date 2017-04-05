function [hyp,p,z] = mengz(varargin)
% This function implements Meng's z-test for correlated correlations (Meng,
% Rubin, & Rosenthal (1992), Comparing Correlated Correlation Coefficients,
% Psych Bulletin 111(1), 172-175.)
%
% mengz(r1, r2, r12, n) compares two correlations r1 and r2:
% r1: correlation between X and Y
% r2: correlation between X and Z
% rx: correlation between Y and Z
% n: number of observations used to compute correlations
%
% mengz(R, k, n) tests the heterogeneity of a correlation matrix with
% respect to correlating with the variable indicated by index k. This test
% is a chi-squared test, so output argument z is actually a chi-squared
% statistic.
%
% mengz(R, k, n, lambda) tests the contrast indicated by vector lambda.
%
% h: hypothesis outcome (1 - reject null hypothesis of equal correlations
% under alpha level of 0.05 (one-tailed))
% p: chance of falsely rejecting null hypothesis
% z: computed z- or chi-square score of Meng's test
%
% Copyright (C) 2012 Eelke Spaak, Donders Institute for Brain, Cognition,
% and Behaviour, Nijmegen, The Netherlands

if nargin == 4 && isscalar(varargin{1})
  % compare two correlations
  
  r1 = varargin{1};
  r2 = varargin{2};
  rx = varargin{3};
  n = varargin{4};
  
  % compute terms
  rsqmean = (r1*r1 + r2*r2)/2;
  f = (1-rx) / (2*(1 - rsqmean));
  if f > 1
    f = 1;
  end
  h = (1-f*rsqmean) / (1 - rsqmean);

  % Fisher transform the two correlations
  z1 = atanh(r1);
  z2 = atanh(r2);

  % compute z
  z = (z1 - z2) * sqrt( (n-3) / (2*(1-rx)*h) );

  % perform one-tailed test
  p = 1-normcdf(z, 0, 1);
  hyp = (p < 0.05);
  
elseif nargin == 3 && isscalar(varargin{2})
  % test heterogeneity of correlation matrix
  
  R = varargin{1};
  k = varargin{2};
  n = varargin{3};
  
  % extract 'predicted' variable from matrix
  X = R(k,:);
  X(k) = [];
  R(k,:) = [];
  R(:,k) = [];
  
  % create vector with all correlations between the predictor variables
  % occurring exactly once
  predr = R(logical(tril(ones(size(R)), -1)));
  
  % compute terms
  rsqmean = mean(X.^2);
  rx = median(predr);
  f = (1-rx) / (2*(1 - rsqmean));
  if f > 1
    f = 1;
  end
  h = (1-f*rsqmean) / (1 - rsqmean);
  
  % Fisher transform
  zX = atanh(X);
  
  % compute chi-square
  z = ( (n-3) * sum( (zX - repmat(mean(zX), size(zX))).^2 ) ) / ( (1-rx)*h );
  
  % perform test
  p = 1-chi2cdf(z, numel(X)-1);
  hyp = (p < 0.05);
  
elseif nargin == 4 && ~isscalar(varargin{1})
  % compute contrast
  
  R = varargin{1};
  k = varargin{2};
  n = varargin{3};
  lambda = varargin{4};
  
  % make sure lambda is a column vector
  if size(lambda,2) > size(lambda,1)
    lambda = lambda';
  end
  
  % use the test of heterogeneity for easy contrast testing
  [~,~,chiSq] = mengz(R,k,n);
  
  % extract 'predicted' variable from matrix and Fisher-transform
  zX = atanh(R(k,:))';
  zX(k) = [];
  
  % compute z-score of contrast
  z = corr(zX, lambda) * sqrt(chiSq);
  
  % perform one-tailed test
  p = 1-normcdf(z, 0, 1);
  hyp = (p < 0.05);
  
else
  error('invalid calling syntax, see HELP MENGZ for more information');
end

end