function [report] = fit_TOJ(TOJdata, LamBounds, TauBounds, DeltaBounds, ...
      LamT_Start, LamR_Start, Tau_Start, Delta_Start, Err_Start, Bias_Start, ...
      Model, Plot, Disp)
%
% Fit model to data from a TOJ task
%
% Usage: o = fit_TOJ(SJ3data, LamBounds, TauBounds, DeltaBounds, ...
%                    LamT_Start, LamR_Start, Tau_Start, ...
%                    Delta_Start, Err_Start, Bias_Start, Model, Plot, Disp)
% Arguments:
%
% - TOJdata:     3-by-N_TOJ array with data from the TOJ task. The first row
%                contains the N_SJ2 stimulus levels; the second row contains the
%                count of TF responses at each level; the third row contains the
%                count of RF responses at each level
%
% - LamBounds:   2-element vector with the lower and upper bounds on lambda_r and
%                lambda_t
%
% - TauBounds:   2-element vector with the lower and upper bounds on tau
%
% - DeltaBounds: 2-element vector with the lower and upper bounds on delta
%
% - LamT_Start:  Starting value(s) for lambda_t. Scalar or vector of positive reals
%
% - LamR_Start:  Starting value(s) for lambda_r. Scalar or vector of positive reals
%
% - Tau_Start:   Starting value(s) for tau. Scalar or vector of reals
%
% - Delta_Start: Starting value(s) for delta. Scalar or vector of positive reals
%
% - Err_Start:   Starting value(s) for epsilon. Scalar or vector of reals in [0, 1]
%
% - Bias_Start:  Starting value(s) for xi and kappa. Scalar or vector of reals in [0, 1]
%
% - Model:       Choice of model to be fitted. An integer scalar in the range between 0
%                and 3 (implying different assumptions about error parameters; see the
%                table underneath) or the string 'best' (case insensitive) to find the
%                best-fitting model by the BIC criterion
%
%                                ___TOJ____
%                    Model #     e_TF  e_RF
%                    -------     ----------
%                          0      x     x
%                          1      o     o
%                          2      o     x
%                          3      x     o
%                    ----------------------------
%                                 x : parameter is not included in the model
%                                 o : parameter is included in the model
%
% - Plot:        Option to plot data and fitted functions. Logical scalar
%
% - Disp:        Option to issue warnings or display progress information. Logical scalar
%
%
% Output: The output is a structure including parameter estimates and goodness-of-fit
%         measures and p-values
%
% Further details in:  Alcalá-Quintana, R. & García-Pérez, M.A. (2013). Fitting model-based
%                         psychometric functions to simultaneity and temporal-order judgment
%                         data: MATLAB and R routines. Behavior Research Methods, in press.
%                         http://dx.doi.org/10.3758/s13428-013-0325-2

% ------------------Check conditions---------------------------------------
if ndims(TOJdata)>2 || size(TOJdata,1)~=3 || ~isnumeric(TOJdata) || ~isreal(TOJdata)
  error('Invalid TOJdata (must be a 2D numeric matrix with 3 rows)')
elseif ~isnumeric(LamBounds) || ~isreal(LamBounds) || any(LamBounds<=0)
  error('Invalid LamBounds (must be positive reals)')
elseif numel(LamBounds)~=2
  error('Invalid size of LamBounds (must have two components only)')
elseif LamBounds(1)==LamBounds(2)
  error('Invalid contents of LamBounds (components must have different values)')
elseif ~isnumeric(TauBounds) || ~isreal(TauBounds)
  error('Invalid TauBounds (must be reals)')
elseif numel(TauBounds)~=2
  error('Invalid size of TauBounds (must have two components only)')
elseif TauBounds(1)==TauBounds(2)
  error('Invalid contents of TauBounds (components must have different values)')
elseif ~isnumeric(DeltaBounds) || ~isreal(DeltaBounds) || any(DeltaBounds<0)
  error('Invalid DeltaBounds (must be positive reals)')
elseif numel(DeltaBounds)~=2
  error('Invalid size of DeltaBounds (must have two components only)')
elseif DeltaBounds(1)==DeltaBounds(2)
  error('Invalid contents of DeltaBounds (components must have different values)')
elseif ~isnumeric(LamT_Start) || ~isreal(LamT_Start)
  error('Invalid LamTStart (must be real)')
elseif ~isnumeric(LamR_Start) || ~isreal(LamR_Start)
  error('Invalid LamRStart (must be real)')
elseif ~isnumeric(Tau_Start) || ~isreal(Tau_Start)
  error('Invalid TauStart (must be real)')
elseif ~isnumeric(Delta_Start) || ~isreal(Delta_Start)
  error('Invalid DeltaStart (must be real)')
elseif ~isnumeric(Err_Start) || ~isreal(Err_Start)
  error('Invalid ErrStart (must be real)')
elseif ~isnumeric(Bias_Start) || ~isreal(Bias_Start)
  error('Invalid BiasStart (must be real)')
elseif any(any(TOJdata(2:3,:)<0))
  error('Negative counts in rows 2-3 of TOJdata')
elseif ~strcmpi(Model,'best') && (~isscalar(Model) || ~isnumeric(Model))
  error('Invalid value for Model (must be the string ''best'' or an integer scalar between 0 and 3)')
elseif isnumeric(Model) && (~isreal(Model) || mod(Model,1)~=0 || Model<0 || Model>3)
  error('Invalid value for Model (must be the string ''best'' or an integer scalar between 0 and 3)')
elseif ~islogical(Plot) || ~isscalar(Plot)
  error('Invalid value for Plot (must be a logical scalar)')
elseif ~islogical(Disp) || ~isscalar(Disp)
  error('Invalid value for Disp (must be a logical scalar)')
end
% ------------------Main---------------------------------------------------
s1 = warning('off','MATLAB:log:logOfZero');
s2 = warning('off','MATLAB:divideByZero');
LamBounds = sort(reshape(LamBounds,1,2));
TauBounds = sort(reshape(TauBounds,1,2));
DeltaBounds = sort(reshape(DeltaBounds,1,2));
LamT_Start = reshape(LamT_Start,1,numel(LamT_Start));
LamR_Start = reshape(LamR_Start,1,numel(LamR_Start));
Tau_Start = reshape(Tau_Start,1,numel(Tau_Start));
Delta_Start = reshape(Delta_Start,1,numel(Delta_Start));
Err_Start = reshape(Err_Start,1,numel(Err_Start));
Bias_Start = reshape(Bias_Start,1,numel(Bias_Start));
LevTOJ = TOJdata(1,:); TOJ_TF = TOJdata(2,:); TOJ_RF = TOJdata(3,:);
NforBIC = sum(TOJ_TF+TOJ_RF>0);
options = optimset('MaxIter',500, 'MaxFunEvals',2000, ...
                   'TolCon',1e-4, 'TolFun',1e-4, 'TolX',1e-4, ...
                   'LargeScale','off', 'Display','off');
rel = version('-release'); rel = sscanf(rel,'%i4');
if rel==2008 || rel==2009, options = optimset(options, 'Algorithm', 'interior-point'); end
if rel>=2010, options = optimset(options, 'Algorithm', 'sqp'); end
if Model==0, Err_Start=0; end
Models = Model; if strcmpi(Model,'best'), Models = 0:3; end; NumModels=numel(Models);
NumIter=length(LamT_Start)*length(LamR_Start)*length(Tau_Start)*length(Delta_Start)* ...
        length(Err_Start)*length(Bias_Start);
ChosenModel = NaN;
BICmin = realmax;
for ModelLoop = Models
  strm = ['Tried ',num2str(ModelLoop+1),' of ',num2str(NumModels),' models; '];
% proceed through initial values
  niter = 0;
  iv_pars = NaN(NumIter,21);
  iv_value = NaN(NumIter,1);
  iv_flag = NaN(NumIter,1);
  iv_iterations = NaN(NumIter,1);
  iv_funcCount = NaN(NumIter,1);
  OnBound = NaN(NumIter,1);
  for LamTInit = LamT_Start
  for LamRInit = LamR_Start
  for TauInit = Tau_Start
  for DeltaInit = Delta_Start
  for ErrInit = Err_Start
  for BiasInit = Bias_Start
    niter = niter+1;
    str = ['Tried ',num2str(niter),'/',num2str(NumIter),' sets of initial values; '];
    if ModelLoop==0
      Initial = [LamTInit LamRInit TauInit DeltaInit BiasInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) 0];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) 1];
    elseif ModelLoop==1
      Initial = [LamTInit LamRInit TauInit DeltaInit BiasInit ErrInit ErrInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) 0 0 0];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) 1 1 1];
    elseif ModelLoop==2 || ModelLoop==3
      Initial = [LamTInit LamRInit TauInit DeltaInit BiasInit ErrInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) 0 0];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) 1 1];
    end
    npars = length(Initial);
    % estimate parameters
    [xmin, value, flag, output] = fmincon(@LogLikelihood, Initial, ...
          [], [], [], [], lb, ub, [], options, ...
          TOJ_TF, TOJ_RF, LevTOJ, ModelLoop);
    iv_pars(niter,1:npars) = xmin;
    iv_value(niter) = value;
    iv_flag(niter) = flag;
    iv_iterations(niter) = output.iterations;
    iv_funcCount(niter) = output.funcCount;
    OnBound(niter) = any(xmin(1:2)==LamBounds(2));
    if NumModels==1 && Disp, disp([str,'Obj.func. = ',num2str(value)]), end
  end
  end
  end
  end
  end
  end
  [v1 best] = min(iv_value);
  if OnBound(best) && any(~OnBound)
    iv_value = iv_value.*(OnBound+1);
    [v2 best2] = min(iv_value);
    if v2<=v1*1.005, best = best2; end
  end
  temp_pars = iv_pars(best,1:npars);
  temp_flag = iv_flag(best);
  temp_iterations = iv_iterations(best);
  temp_funcCount = iv_funcCount(best);
  ObjectiveFunctionMin = iv_value(best);
  BIC = 2*ObjectiveFunctionMin + npars*log(NforBIC);
  if BIC<BICmin
    BICmin=BIC;
    ChosenModel = ModelLoop;
    Final_flag = temp_flag;
    Final_npars = npars;
    Final_pars = temp_pars;
    Final_iterations = temp_iterations;
    Final_funcCount = temp_funcCount;
    if NumModels>1 && Disp
      disp([strm,' BIC = ',num2str(BIC),' ; best thus far ...'])
    end
  elseif NumModels>1 && Disp
    disp([strm,' BIC = ',num2str(BIC)])
  end
end
npars = Final_npars; pars = Final_pars; BIC = BICmin;
% extract estimated parameters
[lam_T lam_R tau delta_TOJ xi epskap] = GetParVec(pars, ChosenModel);
if Disp
  msg1 = 'Consider re-running with a smaller lower bound in LamBounds';
  msg2 = 'Potentially non-informative data';
  if lam_T==LamBounds(1), disp(['Warning: Estimated lambda_T at the lower bound. ' msg1]),end
  if lam_R==LamBounds(1), disp(['Warning: Estimated lambda_R at the lower bound. ' msg1]),end
  if lam_T==LamBounds(2), disp(['Warning: Estimated lambda_T at the upper bound. ' msg2]),end
  if lam_R==LamBounds(2), disp(['Warning: Estimated lambda_R at the upper bound. ' msg2]),end
end
% -----------------Goodness of fit test------------------------------------
% observed frequencies
O = [TOJ_TF TOJ_RF];
% number of observations per delay
n_TOJ = TOJ_TF + TOJ_RF;
% theoretical probabilities
P_TOJ_TF = psi_TOJ_TF(LevTOJ, lam_T, lam_R, tau, delta_TOJ, xi, epskap);
P_TOJ_RF = 1 - P_TOJ_TF;
% expectations
E = [P_TOJ_TF.*n_TOJ P_TOJ_RF.*n_TOJ];
% degrees of freedom
dof = sum(n_TOJ>0) - npars;
% index of fit and significance
EE = E(E>0); OO = O(E>0);
G2 = 2*nansum(OO.*log(OO./EE));
X2 = sum((OO-EE).^2./EE);
pG2 = 'Insufficient degrees of freedom';
pX2 = 'Insufficient degrees of freedom';
if dof > 0
  pG2 = 1 - chi2cdf(G2, dof);
  pX2 = 1 - chi2cdf(X2, dof);
end
if Plot, PlotIt(pars, LevTOJ, ChosenModel, TOJ_TF, TOJ_RF); end
% rearrange parameters to return
params = [pars(1:4) epskap(1) xi epskap(2)];
if ChosenModel==0, params(5)=NaN; params(7)=NaN;
elseif ChosenModel==2, params(7)=NaN;
elseif ChosenModel==3, params(5)=NaN;
end
% create report
report = struct('Problem', ['Separate fit of TOJ data with ' output.algorithm],...
                'FlagFrom_fmincon', Final_flag, ...
                'NumIterations', Final_iterations, ...
                'NumFuncCount', Final_funcCount, ...
                'UserModel', Model, ...
                'FittedModel', ChosenModel, ...
                'NumFreeParameters', npars, ...
                'NumCells', 2*sum(n_TOJ>0), ...
                'NumExpBelow5', sum(E<5), ...
                'NumExpBelow5_ObsAbove0', sum((E<5).*(O>0)), ...
                'NumExpBelow1', sum(E<1), ...
                'NumExpBelow1_ObsAbove0', sum((E<1).*(O>0)), ...
                'DegreesOfFreedom', dof, ...
                'ChiSquareStatistic', X2, ...
                'ChiSquareP_value', pX2, ...
                'LikelihoodRatioStatistic', G2, ...
                'LikelihoodRatioP_value', pG2, ...
                'BIC', BIC, ...
                'LamBounds', LamBounds, ...
                'TauBounds', TauBounds, ...
                'DeltaBounds', DeltaBounds, ...
                'TOJ_LamT', params(1), ...
                'TOJ_LamR', params(2), ...
                'TOJ_Tau', params(3), ...
                'TOJ_Delta', params(4), ...
                'TOJ_epsilonTF', params(5), ...
                'TOJ_Xi', params(6), ...
                'TOJ_epsilonRF', params(7));
warning(s1); warning(s2);
% -----------------Calculate Likelihood------------------------------------
function L = LogLikelihood(x, TOJ_TF, TOJ_RF, LevTOJ, Model)
  [lam_T lam_R tau delta_TOJ xi epskap] = GetParVec(x, Model);
  penalty = 1.0e20;
  if lam_T<=0 || lam_R<=0 || delta_TOJ<0 || any([xi epskap]<0) || any([xi epskap]>1)
    L = penalty;
  else
    P_TOJ_TF = psi_TOJ_TF(LevTOJ, lam_T, lam_R, tau, delta_TOJ, xi, epskap);
    P_TOJ_RF = 1 - P_TOJ_TF;
    L = -nansum([TOJ_TF.*log(P_TOJ_TF) TOJ_RF.*log(P_TOJ_RF)]);
    if ~isreal(L) || isinf(L), L = penalty; end
  end
% -----------------Psychometric Functions----------------------------------
function p = psi_TOJ_TF(x, lam_T, lam_R, tau, delta, xi, epskap)
  lTF = epskap(1); lRF = epskap(2);
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = (1-lTF)*pTF + (1-xi)*pS + lRF*pRF;
% -----------------Get parameters out of vector of estimates---------------
function [lam_T lam_R tau delta_TOJ xi epskap] = GetParVec(p, Model)
lam_T = p(1);  lam_R = p(2);  tau = p(3); delta_TOJ = p(4); xi = p(5);
if Model==0, epskap = zeros(1,2);
elseif Model==1, epskap = p(6:7);
elseif Model==2, epskap = [p(6) 0];
elseif Model==3, epskap = [0 p(6)];
end
%------------------Distribution of the arrival-time-difference-------------
function d = difcdf(x,shift,lam1,lam2)
  y = x-shift;
  left = min(lam1.*exp(lam2.*y)./(lam1+lam2),1);
  right = max(1-(lam2.*exp(-lam1.*y))./(lam1+lam2),0);
  d = (y<=0).*left + (y>0).*right;
% -----------------Plot Results--------------------------------------------
function PlotIt(pars, LevTOJ, Model, TOJ_TF, TOJ_RF)
[lam_T lam_R tau delta_TOJ xi epskap] = GetParVec(pars, Model);
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/10 scrsz(4)/3 4*scrsz(3)/15 scrsz(4)/4])
lw = 3; % line width
ms = 8; % marker size
fsa = 12; % font size for axis labels and title
FS = ['\fontsize{',num2str(fsa),'}'];
fsl = 11; % font size for legend
  subplot (1,1,1)
xinf = min(LevTOJ);  xsup = max(LevTOJ);
step = (xsup-xinf)/800;
xval = xinf:step:xsup; % values x-axis
P_TOJ_TF = psi_TOJ_TF(xval, lam_T, lam_R, tau, delta_TOJ, xi, epskap);
P_TOJ_RF = 1 - P_TOJ_TF;
n = TOJ_TF + TOJ_RF;
plot(xval,P_TOJ_RF,'c','LineWidth',lw); hold on
plot(LevTOJ,5./n,'k -.',LevTOJ,1./n,'k :','LineWidth',1)
plot(LevTOJ,TOJ_RF./n,'cO','MarkerFaceColor','c','MarkerSize',ms); hold off
axis([xinf xsup 0 1])
xlabel([FS ' Delay or SOA, \Delta\itt'])
ylabel([FS ' Probability of response'])
title([FS ' TOJ task'])
text(xinf+0.05*(xsup-xinf),0.9,'\Psi^{*}_{TOJ-RF}',...
     'Color','c','FontWeight','demi','FontSize',fsl)
