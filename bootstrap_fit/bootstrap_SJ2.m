function [out] = bootstrap_SJ2(report, SJ2data, SampleSize, ConfCoef, FixedSeed)
%
% Obtain bootstrap confidence intervals (CIs) and p-values for model fitted
% to data from an SJ2 task through fit_SJ2
%
% Usage: o = bootstrap_SJ2(report, SJ2data, ...
%                          SampleSize, ConfCoef, FixedSeed)
% Arguments:
%
% - report:     Output structure returned by fit_SJ2
%
% - SJ2data:    Data matrix given to fit_SJ2, which produced the output structure
%               entered as the first argument
%
% - SampleSize: Requested umber of bootstrap samples (integer scalar greater than 1)
%
% - ConfCoef:   Confidence coefficient for CIs (real scalar between 80 and 99.9)
%
% - FixeSeed:   Option to fix the seed of the random number generator (logical scalar)
%
%
% Output: The output is a structure including CIs and bootstrap p-values
%
% Further details in:  Alcalá-Quintana, R. & García-Pérez, M.A. (2013). Fitting model-based
%                         psychometric functions to simultaneity and temporal-order judgment
%                         data: MATLAB and R routines. Behavior Research Methods, in press.
%                         http://dx.doi.org/10.3758/s13428-013-0325-2

% ------------------Check conditions---------------------------------------
if ndims(SJ2data)>2 || size(SJ2data,1)~=3 || ~isnumeric(SJ2data) || ~isreal(SJ2data)
  error('Invalid SJ2data (must be a 2D numeric matrix with 3 rows)')
elseif any(any(SJ2data(2:3,:)<0))
  error('Negative counts in rows 2-3 of SJ2data')
elseif ~isstruct(report)
  error('First argument must be an output structure from fit_SJ2')
elseif ~strncmpi(report.Problem,'Separate fit of SJ2',19)
  error('Input structure does not seem to come from fit_SJ2')
elseif ~islogical(FixedSeed) || ~isscalar(FixedSeed)
  error('Invalid value for FixedSeed (must be a logical scalar)')
elseif ~isnumeric(SampleSize) || ~isreal(SampleSize) || mod(SampleSize,1)~=0 || SampleSize<=1 || numel(SampleSize)~=1
  error('Invalid SampleSize (must be an integer greater than 1)')
elseif ~isnumeric(ConfCoef) || ~isreal(ConfCoef) || ConfCoef<80 || ConfCoef>99.9 || numel(ConfCoef)~=1
  error('Invalid ConfCoef (must be a real scalar between 80 and 99.9)')
end
% ------------------Main---------------------------------------------------
if FixedSeed
  seed=931316785;
  rand('state',seed);
  randn('state',seed);
end
estimates = NaN(SampleSize,7);
stat = NaN(SampleSize,2);
s1 = warning('off','MATLAB:log:logOfZero');
s2 = warning('off','MATLAB:divideByZero');
pars = [report.SJ2_LamT report.SJ2_LamR report.SJ2_Tau ...
        report.SJ2_Delta report.SJ2_epsilonTF report.SJ2_epsilonS report.SJ2_epsilonRF];
PE = pars';
pars(isnan(pars)) = 0;
Model = report.FittedModel;
LamBounds = report.LamBounds;
TauBounds = report.TauBounds;
DeltaBounds = report.DeltaBounds;
LamT_Start = pars(1);
LamR_Start = pars(2);
Tau_Start = pars(3);
Delta_Start = pars(4);
Err_Start = 0.05;
LevSJ2 = SJ2data(1,:); Trials_SJ2 = SJ2data(2,:)+SJ2data(3,:);
SJ2dat(1,:) = LevSJ2;
epskap = pars(5:7);
howmany = SampleSize/10;
for i=1:SampleSize
  if mod(i,howmany)==0, disp([num2str(i),' of ',num2str(SampleSize),' samples drawn ...']), end
  % simulate SJ2
  p1 = psi_SJ2_S(LevSJ2,pars(1),pars(2),pars(3),pars(4),epskap);
  sim = binornd(Trials_SJ2,p1);
  SJ2dat(3,:) = sim;
  SJ2dat(2,:) = Trials_SJ2-sim;
  t = fit_SJ2(SJ2dat, LamBounds, TauBounds, DeltaBounds, ...
              LamT_Start, LamR_Start, Tau_Start, Delta_Start, Err_Start, ...
              Model, false, false);
  estimates(i,:) = [t.SJ2_LamT t.SJ2_LamR t.SJ2_Tau ...
                    t.SJ2_Delta t.SJ2_epsilonTF t.SJ2_epsilonS t.SJ2_epsilonRF];
  stat(i,1) = t.ChiSquareStatistic;
  stat(i,2) = t.LikelihoodRatioStatistic;
end
alpha2 = (1-ConfCoef/100)/2;
CI = [PE prctile(estimates,100*[alpha2 1-alpha2])'];
p_value_X2 = sum(stat(:,1)>report.ChiSquareStatistic)/SampleSize;
p_value_G2 = sum(stat(:,2)>report.LikelihoodRatioStatistic)/SampleSize;
% create report
out = struct('Problem', 'Bootstrap results for fit_SJ2',...
             'Model', Model, ...
             'SampleSize', SampleSize, ...
             'ConfidenceLevel', ConfCoef, ...
             'PE_CI_SJ2_LamT', CI(1,:), ...
             'PE_CI_SJ2_LamR', CI(2,:), ...
             'PE_CI_SJ2_Tau', CI(3,:), ...
             'PE_CI_SJ2_Delta', CI(4,:), ...
             'PE_CI_SJ2_epsilonTF', CI(5,:), ...
             'PE_CI_SJ2_epsilonS', CI(6,:), ...
             'PE_CI_SJ2_epsilonRF', CI(7,:), ...
             'pvalue_X2', p_value_X2, ...
             'pvalue_G2', p_value_G2);
warning(s1); warning(s2);
% -----------------Psychometric Functions----------------------------------
function p = psi_SJ2_S(x, lam_T, lam_R, tau, delta, epskap)
  lTF = epskap(1); lS = epskap(2); lRF = epskap(3);
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = lTF*pTF + (1-lS)*pS + lRF*pRF;
%------------------Distribution of the arrival-time difference-------------
function d = difcdf(x,shift,lam1,lam2)
  y = x-shift;
  left = min(lam1.*exp(lam2.*y)./(lam1+lam2),1);
  right = max(1-(lam2.*exp(-lam1.*y))./(lam1+lam2),0);
  d = (y<=0).*left + (y>0).*right;
