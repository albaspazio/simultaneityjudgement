function [out] = bootstrap_TOJ(report, TOJdata, SampleSize, ConfCoef, FixedSeed)
%
% Obtain bootstrap confidence intervals (CIs) and p-values for model fitted
% to data from a TOJ task through fit_TOJ
%
% Usage: o = bootstrap_TOJ(report, TOJdata, ...
%                          SampleSize, ConfCoef, FixedSeed)
% Arguments:
%
% - report:     Output structure returned by fit_TOJ
%
% - TOJdata:    Data matrix given to fit_TOJ, which produced the output structure
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
if ndims(TOJdata)>2 || size(TOJdata,1)~=3 || ~isnumeric(TOJdata) || ~isreal(TOJdata)
  error('Invalid TOJdata (must be a 2D numeric matrix with 3 rows)')
elseif any(any(TOJdata(2:3,:)<0))
  error('Negative counts in rows 2-3 of TOJdata')
elseif ~isstruct(report)
  error('First argument must be an output structure from fit_TOJ')
elseif ~strncmpi(report.Problem,'Separate fit of TOJ',19)
  error('Input structure does not seem to come from fit_TOJ') 
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
pars = [report.TOJ_LamT report.TOJ_LamR report.TOJ_Tau ...
        report.TOJ_Delta report.TOJ_epsilonTF report.TOJ_Xi report.TOJ_epsilonRF];
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
Bias_Start = 0.5;
LevTOJ = TOJdata(1,:); Trials_TOJ = TOJdata(2,:)+TOJdata(3,:);
TOJdat(1,:) = LevTOJ;
epskap = [pars(5) pars(7)];
howmany = SampleSize/10;
for i=1:SampleSize
  if mod(i,howmany)==0, disp([num2str(i),' of ',num2str(SampleSize),' samples drawn ...']), end
  p1 = psi_TOJ_TF(LevTOJ,pars(1),pars(2),pars(3),pars(4),pars(6),epskap);
  sim = binornd(Trials_TOJ,p1);
  TOJdat(2,:) = sim;
  TOJdat(3,:) = Trials_TOJ-sim;
  t = fit_TOJ(TOJdat, LamBounds, TauBounds, DeltaBounds, ...
              LamT_Start, LamR_Start, Tau_Start, Delta_Start, Err_Start, Bias_Start, ...
              Model, false, false);
  estimates(i,:) = [t.TOJ_LamT t.TOJ_LamR t.TOJ_Tau ...
                    t.TOJ_Delta t.TOJ_epsilonTF t.TOJ_Xi t.TOJ_epsilonRF];
  stat(i,1) = t.ChiSquareStatistic;
  stat(i,2) = t.LikelihoodRatioStatistic;
end
alpha2 = (1-ConfCoef/100)/2;
CI = [PE prctile(estimates,100*[alpha2 1-alpha2])'];
p_value_X2 = sum(stat(:,1)>report.ChiSquareStatistic)/SampleSize;
p_value_G2 = sum(stat(:,2)>report.LikelihoodRatioStatistic)/SampleSize;
% create report
out = struct('Problem', 'Bootstrap results for fit_TOJ',...
             'Model', Model, ...
             'SampleSize', SampleSize, ...
             'ConfidenceLevel', ConfCoef, ...
             'PE_CI_TOJ_LamT', CI(1,:), ...
             'PE_CI_TOJ_LamR', CI(2,:), ...
             'PE_CI_TOJ_Tau', CI(3,:), ...
             'PE_CI_TOJ_Delta', CI(4,:), ...
             'PE_CI_TOJ_epsilonTF', CI(5,:), ...
             'PE_CI_TOJ_Xi', CI(6,:), ...
             'PE_CI_TOJ_epsilonRF', CI(7,:), ...
             'pvalue_X2', p_value_X2, ...
             'pvalue_G2', p_value_G2);
warning(s1); warning(s2);
% -----------------Psychometric Functions----------------------------------
function p = psi_TOJ_TF(x, lam_T, lam_R, tau, delta, xi, epskap)
  lTF = epskap(1); lRF = epskap(2);
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = (1-lTF)*pTF + (1-xi)*pS + lRF*pRF;
%------------------Distribution of the arrival-time difference-------------
function d = difcdf(x,shift,lam1,lam2)
  y = x-shift;
  left = min(lam1.*exp(lam2.*y)./(lam1+lam2),1);
  right = max(1-(lam2.*exp(-lam1.*y))./(lam1+lam2),0);
  d = (y<=0).*left + (y>0).*right;
