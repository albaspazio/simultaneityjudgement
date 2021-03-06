function [out] = bootstrap_SJ3_TOJ(report, SJ3data, TOJdata, ...
                                   SampleSize, ConfCoef, FixedSeed)
%
% Obtain bootstrap confidence intervals (CIs) and p-values for model fitted jointly
% to data from SJ3 and TOJ tasks through fit_SJ3_TOJ
%
% Usage: o = bootstrap_SJ3_TOJ(report, SJ3data, TOJdata, ...
%                              SampleSize, ConfCoef, FixedSeed)
% Arguments:
%
% - report:     Output structure returned by fit_SJ3_TOJ
%
% - SJ3data,
%   TOJdata:    Data matrices given to fit_SJ3_TOJ, which produced the output structure
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
% Further details in:  Alcal?-Quintana, R. & Garc?a-P?rez, M.A. (2013). Fitting model-based
%                         psychometric functions to simultaneity and temporal-order judgment
%                         data: MATLAB and R routines. Behavior Research Methods, in press.
%                         http://dx.doi.org/10.3758/s13428-013-0325-2

% ------------------Check conditions---------------------------------------
if ndims(SJ3data)>2 || size(SJ3data,1)~=4 || ~isnumeric(SJ3data) || ~isreal(SJ3data)
  error('Invalid SJ3data (must be a 2D numeric matrix with 4 rows)')
elseif ndims(TOJdata)>2 || size(TOJdata,1)~=3 || ~isnumeric(TOJdata) || ~isreal(TOJdata)
  error('Invalid TOJdata (must be a 2D numeric matrix with 3 rows)')
elseif any(any(SJ3data(2:4,:)<0))
  error('Negative counts in rows 2-4 of SJ3data')
elseif any(any(TOJdata(2:3,:)<0))
  error('Negative counts in rows 2-3 of TOJdata')
elseif ~isstruct(report)
  error('First argument must be an output structure from fit_SJ3_TOJ')
elseif ~strncmpi(report.Problem,'Joint fit of SJ3 and TOJ',24)
  error('Input structure does not seem to come from fit_SJ3_TOJ')
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
estimates = NaN(SampleSize,17);
stat = NaN(SampleSize,2);
s1 = warning('off','MATLAB:log:logOfZero');
s2 = warning('off','MATLAB:divideByZero');
pars = [report.common_LamT report.common_LamR report.common_Tau ...
        report.SJ3_Delta report.SJ3_epsilonTF report.SJ3_epsilonS report.SJ3_epsilonRF ...
        report.SJ3_kappa_TFintoS report.SJ3_kappa_TFintoRF ...
        report.SJ3_kappa_SintoTF report.SJ3_kappa_SintoRF ...
        report.SJ3_kappa_RFintoTF report.SJ3_kappa_RFintoS ...
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
Delta_Start = mean([pars(4) pars(14)]);
Err_Start = 0.05;
Bias_Start = 0.5;
LevSJ3 = SJ3data(1,:); Trials_SJ3 = SJ3data(2,:)+SJ3data(3,:)+SJ3data(4,:);
LevTOJ = TOJdata(1,:); Trials_TOJ = TOJdata(2,:)+TOJdata(3,:);
SJ3dat(1,:) = LevSJ3; TOJdat(1,:) = LevTOJ;
epskap = [pars(5:13) pars(15) pars(17)];
howmany = SampleSize/10;
for i=1:SampleSize
  if mod(i,howmany)==0, disp([num2str(i),' of ',num2str(SampleSize),' samples drawn ...']), end
  % simulate SJ3
  [p1 p2] = psi_SJ3_TF_S(LevSJ3,pars(1),pars(2),pars(3),pars(4),epskap(1:9));
  p3 = 1-p1-p2;
  sim1 = binornd(Trials_SJ3,p1);
  sim2 = binornd(Trials_SJ3-sim1,p2./(p2+p3));
  SJ3dat(2,:) = sim1;
  SJ3dat(3,:) = sim2;
  SJ3dat(4,:) = Trials_SJ3-sim1-sim2;
  % simulate TOJ
  p1 = psi_TOJ_TF(LevTOJ,pars(1),pars(2),pars(3),pars(14),pars(16),epskap(10:11));
  sim = binornd(Trials_TOJ,p1);
  TOJdat(2,:) = sim;
  TOJdat(3,:) = Trials_TOJ-sim;
  t = fit_SJ3_TOJ(SJ3dat, TOJdat, LamBounds, TauBounds, DeltaBounds, ...
                  LamT_Start, LamR_Start, Tau_Start, Delta_Start, Err_Start, Bias_Start, ...
                  Model, false, false);
  estimates(i,:) = [t.common_LamT t.common_LamR t.common_Tau ...
                    t.SJ3_Delta t.SJ3_epsilonTF t.SJ3_epsilonS t.SJ3_epsilonRF ...
                    t.SJ3_kappa_TFintoS t.SJ3_kappa_TFintoRF ...
                    t.SJ3_kappa_SintoTF t.SJ3_kappa_SintoRF ...
                    t.SJ3_kappa_RFintoTF t.SJ3_kappa_RFintoS ...
                    t.TOJ_Delta t.TOJ_epsilonTF t.TOJ_Xi t.TOJ_epsilonRF];
  stat(i,1) = t.ChiSquareStatistic;
  stat(i,2) = t.LikelihoodRatioStatistic;
end
alpha2 = (1-ConfCoef/100)/2;
CI = [PE prctile(estimates,100*[alpha2 1-alpha2])'];
p_value_X2 = sum(stat(:,1)>report.ChiSquareStatistic)/SampleSize;
p_value_G2 = sum(stat(:,2)>report.LikelihoodRatioStatistic)/SampleSize;
% create report
out = struct('Problem', 'Bootstrap results for fit_SJ3_TOJ',...
             'Model', Model, ...
             'SampleSize', SampleSize, ...
             'ConfidenceLevel', ConfCoef, ...
             'PE_CI_common_LamT', CI(1,:), ...
             'PE_CI_common_LamR', CI(2,:), ...
             'PE_CI_common_Tau', CI(3,:), ...
             'PE_CI_SJ3_Delta', CI(4,:), ...
             'PE_CI_SJ3_epsilonTF', CI(5,:), ...
             'PE_CI_SJ3_epsilonS', CI(6,:), ...
             'PE_CI_SJ3_epsilonRF', CI(7,:), ...
             'PE_CI_SJ3_kappa_TFintoS', CI(8,:), ...
             'PE_CI_SJ3_kappa_TFintoRF', CI(9,:), ...
             'PE_CI_SJ3_kappa_SintoTF', CI(10,:), ...
             'PE_CI_SJ3_kappa_SintoRF', CI(11,:), ...
             'PE_CI_SJ3_kappa_RFintoTF', CI(12,:), ...
             'PE_CI_SJ3_kappa_RFintoS', CI(13,:), ...
             'PE_CI_TOJ_Delta', CI(14,:), ...
             'PE_CI_TOJ_epsilonTF', CI(15,:), ...
             'PE_CI_TOJ_Xi', CI(16,:), ...
             'PE_CI_TOJ_epsilonRF', CI(17,:), ...
             'pvalue_X2', p_value_X2, ...
             'pvalue_G2', p_value_G2);
warning(s1); warning(s2);
% -----------------Psychometric Functions----------------------------------
function [p1 p2] = psi_SJ3_TF_S(x, lam_T, lam_R, tau, delta, epskap)
  lTF = epskap(1); lS = epskap(2); lRF = epskap(3);
  kSTF = epskap(6); kRFTF = epskap(8);
  kTFS = epskap(4); kRFS = epskap(9);
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p1 = (1-lTF)*pTF + lS*kSTF*pS + lRF*kRFTF*pRF;
  p2 = lTF*kTFS*pTF + (1-lS)*pS + lRF*kRFS*pRF;
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
