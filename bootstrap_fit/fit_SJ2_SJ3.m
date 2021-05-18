function [report] = fit_SJ2_SJ3(SJ2data, SJ3data, LamBounds, TauBounds, DeltaBounds, ...
      LamT_Start, LamR_Start, Tau_Start, Delta_Start, Err_Start, Bias_Start, ...
      Model, Plot, Disp)
%
% Fit model jointly to data from SJ2 and SJ3 tasks
%
% Usage: o = fit_SJ2_SJ3(SJ2data, SJ3data, LamBounds, TauBounds, DeltaBounds, ...
%                        LamT_Start, LamR_Start, Tau_Start, ...
%                        Delta_Start, Err_Start, Bias_Start, ...
%                        Model, Plot, Disp)
% Arguments:
%
% - SJ2data:     3-by-N_SJ2 array with data from the SJ2 task. The first row
%                contains the N_SJ2 stimulus levels; the second row contains the
%                count of A responses at each level; the third row contains the
%                count of S responses at each level
%
% - SJ3data:     4-by-N_SJ3 array with data from the SJ3 task. The first row
%                contains the N_SJ3 stimulus levels; the second row contains the
%                count of TF responses at each level; the third row contains the
%                count of S responses at each level; the fourth row contains the
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
% - Model:       Choice of model to be fitted. In the general form, this is a 2-element
%                vector whose components indicate the model chosen for the SJ2 and SJ3
%                tasks. The cases Model = [0 0] and Model = [1 1] can be stated simply
%                as Model = 0 or Model = 1. Valid models (implying different assumptions
%                about error parameters) are given in the table underneath.
%                Alternatively, Model = 'best' (case insensitive) finds the best-fitting
%                model by the BIC criterion
%
%                           ______SJ2_______    ______SJ3_______
%                Model #    e_TF  e_S   e_RF    e_TF  e_S   e_RF
%                -------    ------------------------------------
%                      0     x     x     x       x     x     x  
%                      1     o     o     o       o     o     o  
%                      2     o     o     x       o     o     x  
%                      3     o     x     o       o     x     o  
%                      4     x     o     o       x     o     o  
%                      5     o     x     x       o     x     x  
%                      6     x     o     x       x     o     x 
%                      7     x     x     o       x     x     o 
%                ----------------------------------------------
%                            x : parameter is not included in the model
%                            o : parameter is included in the model
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
if ndims(SJ2data)>2 || size(SJ2data,1)~=3 || ~isnumeric(SJ2data) || ~isreal(SJ2data)
  error('Invalid SJ2data (must be a 2D numeric matrix with 3 rows)')
elseif ndims(SJ3data)>2 || size(SJ3data,1)~=4 || ~isnumeric(SJ3data) || ~isreal(SJ3data)
  error('Invalid SJ3data (must be a 2D numeric matrix with 4 rows)')
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
elseif any(any(SJ2data(2:3,:)<0))
  error('Negative counts in rows 2-3 of SJ2data')
elseif any(any(SJ3data(2:4,:)<0))
  error('Negative counts in rows 2-4 of SJ3data')
elseif ischar(Model) && ~strcmpi(Model,'best')
  error('Invalid string for Model (the only valid string is ''best'')')
elseif ~ischar(Model)
  if isnumeric(Model)
    if numel(Model)==1 && Model~=0 && Model~=1
      error('Invalid value for Model')
    elseif numel(Model)==2 && (any(mod(Model,1))~=0 || any(Model<0) || any(Model>7))
      error('Invalid value for Model')
    elseif numel(Model)~=1 && numel(Model)~=2
      error('Invalid value for Model')
    end
  else
    error('Invalid value for Model (must be numeric or string)')   
  end
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
LevSJ2 = SJ2data(1,:); SJ2_A = SJ2data(2,:); SJ2_S = SJ2data(3,:);
LevSJ3 = SJ3data(1,:); SJ3_TF = SJ3data(2,:); SJ3_S = SJ3data(3,:); SJ3_RF = SJ3data(4,:);
NforBIC = sum(SJ2_A+SJ2_S>0)+2*sum(SJ3_TF+SJ3_S+SJ3_RF>0);
options = optimset('MaxIter',500, 'MaxFunEvals',2000, ...
                   'TolCon',1e-4, 'TolFun',1e-4, 'TolX',1e-4, ...
                   'LargeScale','off', 'Display','off');
rel = version('-release'); rel = sscanf(rel,'%i4');
if rel==2008 || rel==2009, options = optimset(options, 'Algorithm', 'interior-point'); end
if rel>=2010, options = optimset(options, 'Algorithm', 'sqp'); end
if all(Model==0), Err_Start=0; end
if ischar(Model)
  ModelsSJ2 = 0:7; ModelsSJ3 = 0:7;
elseif numel(Model)==1
  ModelsSJ2 = Model; ModelsSJ3 = Model;
else
  ModelsSJ2 = Model(1); ModelsSJ3 = Model(2);
end
NumModels=numel(ModelsSJ2)*numel(ModelsSJ3);
NumIter=length(LamT_Start)*length(LamR_Start)*length(Tau_Start)*length(Delta_Start)* ...
        length(Err_Start)*length(Bias_Start);
ChosenModel = NaN;
BICmin = realmax;
nm = 0;
for ModelLoopSJ2 = ModelsSJ2
for ModelLoopSJ3 = ModelsSJ3
  nm = nm+1;
  strm = ['Tried ',num2str(nm),' of ',num2str(NumModels),' models; '];
  % proceed through initial values
  niter = 0;
  iv_pars = NaN(NumIter,21);
  iv_value = NaN(NumIter,1);
  iv_flag = NaN(NumIter,1);
  iv_iterations = NaN(NumIter,1);
  iv_funcCount = NaN(NumIter,1);
  OnBound = NaN(NumIter,1);
  ObjectiveFunctionMin = realmax;
  for LamTInit = LamT_Start
  for LamRInit = LamR_Start
  for TauInit = Tau_Start
  for DeltaInit = Delta_Start
  for ErrInit = Err_Start
  for BiasInit = Bias_Start
    niter = niter+1;
    str = ['Tried ',num2str(niter),'/',num2str(NumIter),' sets of initial values; '];
    Initial = [LamTInit LamRInit TauInit DeltaInit DeltaInit];
    lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) DeltaBounds(1)];
    ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) DeltaBounds(2)];
    if ModelLoopSJ2==0
    elseif ModelLoopSJ2==1
      Initial = horzcat(Initial,[ErrInit ErrInit ErrInit]);
      lb = horzcat(lb,[0 0 0]); ub = horzcat(ub,[1 1 1]);
    elseif ModelLoopSJ2==2 || ModelLoopSJ2==3 || ModelLoopSJ2==4
      Initial = horzcat(Initial,[ErrInit ErrInit]);
      lb = horzcat(lb,[0 0]); ub = horzcat(ub,[1 1]);
    elseif ModelLoopSJ2==5 || ModelLoopSJ2==6 || ModelLoopSJ2==7
      Initial = horzcat(Initial,ErrInit);
      lb = horzcat(lb,0); ub = horzcat(ub,1);
    end
    if ModelLoopSJ3==0
    elseif ModelLoopSJ3==1
      Initial = horzcat(Initial,[ErrInit ErrInit ErrInit BiasInit BiasInit BiasInit]);
      lb = horzcat(lb,[0 0 0 0 0 0]); ub = horzcat(ub,[1 1 1 1 1 1]);
    elseif ModelLoopSJ3==2 || ModelLoopSJ3==3 || ModelLoopSJ3==4
      Initial = horzcat(Initial,[ErrInit ErrInit BiasInit BiasInit]);
      lb = horzcat(lb,[0 0 0 0]); ub = horzcat(ub,[1 1 1 1]);
    elseif ModelLoopSJ3==5 || ModelLoopSJ3==6 || ModelLoopSJ3==7
      Initial = horzcat(Initial,[ErrInit BiasInit]);
      lb = horzcat(lb,[0 0]); ub = horzcat(ub,[1 1]);
    end
    npars = length(Initial);
    % estimate parameters
    [xmin, value, flag, output] = fmincon(@LogLikelihood, Initial, ...
          [], [], [], [], lb, ub, [], options, ...
          SJ2_S, SJ2_A, SJ3_TF, SJ3_S, SJ3_RF, LevSJ2, LevSJ3, [ModelLoopSJ2 ModelLoopSJ3]);
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
    ChosenModel = [ModelLoopSJ2 ModelLoopSJ3];
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
end
npars = Final_npars; pars = Final_pars; BIC = BICmin;
% extract estimated parameters
[lam_T lam_R tau delta_SJ2 delta_SJ3 epskap] = GetParVec(pars, ChosenModel);
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
O = [SJ2_S SJ2_A SJ3_TF SJ3_S SJ3_RF];
% number of observations per task and delay
n_SJ2 = SJ2_S + SJ2_A;
n_SJ3 = SJ3_TF + SJ3_S + SJ3_RF;
% theoretical probabilities
P_SJ2_S = psi_SJ2_S(LevSJ2, lam_T, lam_R, tau, delta_SJ2, epskap(1:3));
P_SJ2_A = 1 - P_SJ2_S;
[P_SJ3_TF P_SJ3_S] = psi_SJ3_TF_S(LevSJ3, lam_T, lam_R, tau, delta_SJ3, epskap(4:12));
P_SJ3_RF = 1 - P_SJ3_TF - P_SJ3_S;
% expectations
E = [P_SJ2_S.*n_SJ2 P_SJ2_A.*n_SJ2 P_SJ3_TF.*n_SJ3 P_SJ3_S.*n_SJ3 P_SJ3_RF.*n_SJ3];
% degrees of freedom
dof = sum(n_SJ2>0) + 2*sum(n_SJ3>0) - npars;
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
if Plot, PlotIt(pars, LevSJ2, LevSJ3, ChosenModel, SJ2_S, SJ2_A, SJ3_TF, SJ3_S, SJ3_RF); end
% rearrange parameters to return
params = [pars(1:4) epskap(1:3) pars(5) epskap(4:12)];
if ChosenModel(1)==0, params(5:7)=NaN;
elseif ChosenModel(1)==2, params(7)=NaN;
elseif ChosenModel(1)==3, params(6)=NaN;
elseif ChosenModel(1)==4, params(5)=NaN;
elseif ChosenModel(1)==5, params(6:7)=NaN;
elseif ChosenModel(1)==6, params(5)=NaN; params(7)=NaN;
elseif ChosenModel(1)==7, params(5:6)=NaN;
end
if ChosenModel(2)==0, params(9:17)=NaN;
elseif ChosenModel(2)==2, params(11)=NaN; params(16:17)=NaN;
elseif ChosenModel(2)==3, params(10)=NaN; params(14:15)=NaN;
elseif ChosenModel(2)==4, params(9)=NaN; params(12:13)=NaN;
elseif ChosenModel(2)==5, params(10:11)=NaN; params(14:17)=NaN;
elseif ChosenModel(2)==6, params(9)=NaN; params(11)=NaN; params(12:13)=NaN; params(16:17)=NaN;
elseif ChosenModel(2)==7, params(9:10)=NaN; params(12:15)=NaN;
end
% create report
report = struct('Problem', ['Joint fit of SJ2 and SJ3 data with ' output.algorithm],...
                'FlagFrom_fmincon', Final_flag, ...
                'NumIterations', Final_iterations, ...
                'NumFuncCount', Final_funcCount, ...
                'UserModel', Model, ...
                'FittedModel', ChosenModel, ...
                'NumFreeParameters', npars, ...
                'NumCells', 2*sum(n_SJ2>0)+3*sum(n_SJ3>0), ...
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
                'common_LamT', params(1), ...
                'common_LamR', params(2), ...
                'common_Tau', params(3), ...
                'SJ2_Delta', params(4), ...
                'SJ2_epsilonTF', params(5), ...
                'SJ2_epsilonS', params(6), ...
                'SJ2_epsilonRF', params(7), ...
                'SJ3_Delta', params(8), ...
                'SJ3_epsilonTF', params(9), ...
                'SJ3_epsilonS', params(10), ...
                'SJ3_epsilonRF', params(11), ...
                'SJ3_kappa_TFintoS', params(12), ...
                'SJ3_kappa_TFintoRF', params(13), ...
                'SJ3_kappa_SintoTF', params(14), ...
                'SJ3_kappa_SintoRF', params(15), ...
                'SJ3_kappa_RFintoTF', params(16), ...
                'SJ3_kappa_RFintoS', params(17));
warning(s1); warning(s2);
% -----------------Calculate Likelihood------------------------------------
function L = LogLikelihood(x, SJ2_S, SJ2_A, SJ3_TF, SJ3_S, SJ3_RF, LevSJ2, LevSJ3, Model)
  [lam_T lam_R tau delta_SJ2 delta_SJ3 epskap] = GetParVec(x, Model);
  penalty = 1.0e20;
  if lam_T<=0 || lam_R<=0 || delta_SJ2<0 || delta_SJ3<0 || ...
     any(epskap<0) || any(epskap>1)
    L = penalty;
  else
    P_SJ2_S = psi_SJ2_S(LevSJ2, lam_T, lam_R, tau, delta_SJ2, epskap(1:3));
    P_SJ2_A = 1 - P_SJ2_S;
    [P_SJ3_TF P_SJ3_S] = psi_SJ3_TF_S(LevSJ3, lam_T, lam_R, tau, delta_SJ3, epskap(4:12));
    P_SJ3_RF = 1 - P_SJ3_TF - P_SJ3_S;
    L = -nansum([SJ2_S.*log(P_SJ2_S) SJ2_A.*log(P_SJ2_A) ...
                 SJ3_TF.*log(P_SJ3_TF) SJ3_S.*log(P_SJ3_S) SJ3_RF.*log(P_SJ3_RF)]);
    if ~isreal(L) || isinf(L), L = penalty; end
  end
% -----------------Psychometric Functions----------------------------------
function p = psi_SJ2_S(x, lam_T, lam_R, tau, delta, epskap)
  lTF = epskap(1); lS = epskap(2); lRF = epskap(3);
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = lTF*pTF + (1-lS)*pS + lRF*pRF;
function [p1 p2] = psi_SJ3_TF_S(x, lam_T, lam_R, tau, delta, epskap)
  lTF = epskap(1); lS = epskap(2); lRF = epskap(3);
  kSTF = epskap(6); kRFTF = epskap(8);
  kTFS = epskap(4); kRFS = epskap(9);
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p1 = (1-lTF)*pTF + lS*kSTF*pS + lRF*kRFTF*pRF;
  p2 = lTF*kTFS*pTF + (1-lS)*pS + lRF*kRFS*pRF;
% -----------------Get parameters out of vector of estimates---------------
function [lam_T lam_R tau delta_SJ2 delta_SJ3 epskap] = GetParVec(p, Model)
ModelSJ2 = Model(1); ModelSJ3 = Model(2);
lam_T = p(1);  lam_R = p(2);  tau = p(3); delta_SJ2 = p(4);  delta_SJ3 = p(5);
if ModelSJ2==0, epskap = zeros(1,3);
elseif ModelSJ2==1, epskap = p(6:8);
elseif ModelSJ2==2, epskap = [p(6) p(7) 0];
elseif ModelSJ2==3, epskap = [p(6) 0 p(7)];
elseif ModelSJ2==4, epskap = [0 p(6) p(7)];
elseif ModelSJ2==5, epskap = [p(6) 0 0];
elseif ModelSJ2==6, epskap = [0 p(6) 0];
elseif ModelSJ2==7, epskap = [0 0 p(6)];
end
ini=6; if ModelSJ2>=5, ini=7; elseif ModelSJ2>=2, ini=8; elseif ModelSJ2==1, ini=9; end
if ModelSJ3==0, epskap = horzcat(epskap,zeros(1,9));
elseif ModelSJ3==1, epskap = horzcat(epskap,[p(ini:ini+3) 1-p(ini+3) p(ini+4) 1-p(ini+4) p(ini+5) 1-p(ini+5)]);
elseif ModelSJ3==2, epskap = horzcat(epskap,[p(ini:ini+1) 0 p(ini+2) 1-p(ini+2) p(ini+3) 1-p(ini+3) 0 0]);
elseif ModelSJ3==3, epskap = horzcat(epskap,[p(ini) 0 p(ini+1) p(ini+2) 1-p(ini+2) 0 0 p(ini+3) 1-p(ini+3)]);
elseif ModelSJ3==4, epskap = horzcat(epskap,[0 p(ini:ini+1) 0 0 p(ini+2) 1-p(ini+2) p(ini+3) 1-p(ini+3)]);
elseif ModelSJ3==5, epskap = horzcat(epskap,[p(ini) 0 0 p(ini+1) 1-p(ini+1) 0 0 0 0]);
elseif ModelSJ3==6, epskap = horzcat(epskap,[0 p(ini) 0 0 0 p(ini+1) 1-p(ini+1) 0 0]);
elseif ModelSJ3==7, epskap = horzcat(epskap,[0 0 p(ini) 0 0 0 0 p(ini+1) 1-p(ini+1)]);
end
%------------------Distribution of the arrival-time difference-------------
function d = difcdf(x,shift,lam1,lam2)
  y = x-shift;
  left = min(lam1.*exp(lam2.*y)./(lam1+lam2),1);
  right = max(1-(lam2.*exp(-lam1.*y))./(lam1+lam2),0);
  d = (y<=0).*left + (y>0).*right;
% -----------------Plot Results--------------------------------------------
function PlotIt(pars, LevSJ2, LevSJ3, Model, SJ2_S, SJ2_A, SJ3_TF, SJ3_S, SJ3_RF)
[lam_T lam_R tau delta_SJ2 delta_SJ3 epskap] = GetParVec(pars, Model);
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/10 scrsz(4)/3 8*scrsz(3)/15 scrsz(4)/4])
lw = 3; % line width
ms = 8; % marker size
fsa = 12; % font size for axis labels and title
FS = ['\fontsize{',num2str(fsa),'}'];
fsl = 11; % font size for legend
  subplot (1,2,1)
xinf = min(LevSJ2);  xsup = max(LevSJ2);
step = (xsup-xinf)/800;
xval = xinf:step:xsup; % values x-axis
P_SJ2_S = psi_SJ2_S(xval, lam_T, lam_R, tau, delta_SJ2, epskap(1:3));
% P_SJ2_A = 1 - P_SJ2_S; % not used
n = SJ2_S + SJ2_A;
plot(xval,P_SJ2_S,'r','LineWidth',lw); hold on
plot(LevSJ2,5./n,'k -.',LevSJ2,1./n,'k :','LineWidth',1)
plot(LevSJ2,SJ2_S./n,'rO','MarkerFaceColor','r','MarkerSize',ms); hold off
axis([xinf xsup 0 1])
text(xinf+0.05*(xsup-xinf),0.9,'\Psi^{*}_{SJ2-S}',...
     'Color','r','FontWeight','demi','FontSize',fsl)
xlabel([FS ' Delay or SOA, \Delta\itt'])
ylabel([FS ' Probability of response'])
title([FS ' SJ2 task'])
  subplot (1,2,2)
xinf = min(LevSJ3);  xsup = max(LevSJ3);
step = (xsup-xinf)/800;
xval = xinf:step:xsup; % values x-axis
[P_SJ3_TF P_SJ3_S] = psi_SJ3_TF_S(xval, lam_T, lam_R, tau, delta_SJ3, epskap(4:12));
P_SJ3_RF = 1 - P_SJ3_TF - P_SJ3_S;
n = SJ3_TF + SJ3_S + SJ3_RF;
plot(xval,P_SJ3_TF,'k','LineWidth',lw); hold on
plot(xval,P_SJ3_S,'r','LineWidth',lw)
plot(xval,P_SJ3_RF,'c','LineWidth',lw)
plot(LevSJ3,5./n,'k -.',LevSJ3,1./n,'k :','LineWidth',1)
plot(LevSJ3,SJ3_TF./n,'kO','MarkerFaceColor','k','MarkerSize',ms)
plot(LevSJ3,SJ3_S./n,'rO','MarkerFaceColor','r','MarkerSize',ms)
plot(LevSJ3,SJ3_RF./n,'cO','MarkerFaceColor','c','MarkerSize',ms); hold off
axis([xinf xsup 0 1]);
xlabel([FS ' Delay or SOA, \Delta\itt'])
ylabel([FS ' Probability of response'])
title([FS ' SJ3 task'])
text(xinf+0.05*(xsup-xinf),0.55,'\Psi^{*}_{SJ3-TF}', ...
     'Color','k','FontWeight','demi','FontSize',fsl)
text(xinf+0.05*(xsup-xinf),0.45,'\Psi^{*}_{SJ3-S}', ...
     'Color','r','FontWeight','demi','FontSize',fsl)
text(xsup-0.01*(xsup-xinf),0.5,'\Psi^{*}_{SJ3-RF}','HorizontalAlignment','Right', ...
     'Color','c','FontWeight','demi','FontSize',fsl)
