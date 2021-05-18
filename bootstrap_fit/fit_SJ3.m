function [report] = fit_SJ3(SJ3data, LamBounds, TauBounds, DeltaBounds, ...
      LamT_Start, LamR_Start, Tau_Start, Delta_Start, Err_Start, Bias_Start, ...
      Model, Plot, Disp)
%
% Fit model to data from an SJ3 task
%
% Usage: o = fit_SJ3(SJ3data, LamBounds, TauBounds, DeltaBounds, ...
%                    LamT_Start, LamR_Start, Tau_Start, ...
%                    Delta_Start, Err_Start, Bias_Start, Model, Plot, Disp)
% Arguments:
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
% - Model:       Choice of model to be fitted. An integer scalar in the range between 0
%                and 7 (implying different assumptions about error parameters; see the
%                table underneath) or the string 'best' (case insensitive) to find the
%                best-fitting model by the BIC criterion
%
%                                ______SJ3_______
%                    Model #     e_TF  e_S   e_RF
%                    -------     ----------------
%                          0      x     x     x
%                          1      o     o     o
%                          2      o     o     x
%                          3      o     x     o
%                          4      x     o     o
%                          5      o     x     x 
%                          6      x     o     x
%                          7      x     x     o
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
if ndims(SJ3data)>2 || size(SJ3data,1)~=4 || ~isnumeric(SJ3data) || ~isreal(SJ3data)
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
elseif any(any(SJ3data(2:4,:)<0))
  error('Negative counts in rows 2-4 of SJ3data')
elseif ~strcmpi(Model,'best') && (~isscalar(Model) || ~isnumeric(Model))
  error('Invalid value for Model (must be the string ''best'' or an integer scalar between 0 and 7)')
elseif isnumeric(Model) && (~isreal(Model) || mod(Model,1)~=0 || Model<0 || Model>7)
  error('Invalid value for Model (must be the string ''best'' or an integer scalar between 0 and 7)')
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
LevSJ3 = SJ3data(1,:); SJ3_TF = SJ3data(2,:); SJ3_S = SJ3data(3,:); SJ3_RF = SJ3data(4,:);
NforBIC = 2*sum(SJ3_TF+SJ3_S+SJ3_RF>0);
options = optimset('MaxIter',500, 'MaxFunEvals',2000, ...
                   'TolCon',1e-4, 'TolFun',1e-4, 'TolX',1e-4, ...
                   'LargeScale','off', 'Display','off');
rel = version('-release'); rel = sscanf(rel,'%i4');
if rel==2008 || rel==2009, options = optimset(options, 'Algorithm', 'interior-point'); end
if rel>=2010, options = optimset(options, 'Algorithm', 'sqp'); end
if Model==0, Err_Start=0; Bias_Start=0; end
Models = Model; if strcmpi(Model,'best'), Models = 0:7; end; NumModels=numel(Models);
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
      Initial = [LamTInit LamRInit TauInit DeltaInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1)];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2)];
    elseif ModelLoop==1
      Initial = [LamTInit LamRInit TauInit DeltaInit ...
                 ErrInit ErrInit ErrInit BiasInit BiasInit BiasInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) 0 0 0 0 0 0];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) 1 1 1 1 1 1];
    elseif ModelLoop==2 || ModelLoop==3 || ModelLoop==4
      Initial = [LamTInit LamRInit TauInit DeltaInit ...
                 ErrInit ErrInit BiasInit BiasInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) 0 0 0 0];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) 1 1 1 1];
    elseif ModelLoop==5 || ModelLoop==6 || ModelLoop==7
      Initial = [LamTInit LamRInit TauInit DeltaInit ...
                 ErrInit BiasInit];
      lb = [LamBounds(1) LamBounds(1) TauBounds(1) DeltaBounds(1) 0 0];
      ub = [LamBounds(2) LamBounds(2) TauBounds(2) DeltaBounds(2) 1 1];
    end
    npars = length(Initial);
    % estimate parameters
    [xmin, value, flag, output] = fmincon(@LogLikelihood, Initial, ...
          [], [], [], [], lb, ub, [], options, ...
          SJ3_TF, SJ3_S, SJ3_RF, LevSJ3, ModelLoop);
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
[lam_T lam_R tau delta_SJ3 epskap] = GetParVec(pars, ChosenModel);
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
O = [SJ3_TF SJ3_S SJ3_RF];
% number of observations per delay
n_SJ3 = SJ3_TF + SJ3_S + SJ3_RF;
% theoretical probabilities
[P_SJ3_TF P_SJ3_S] = psi_SJ3_TF_S(LevSJ3, lam_T, lam_R, tau, delta_SJ3, epskap);
P_SJ3_RF = 1 - P_SJ3_TF - P_SJ3_S;
% expectations
E = [P_SJ3_TF.*n_SJ3 P_SJ3_S.*n_SJ3 P_SJ3_RF.*n_SJ3];
% degrees of freedom
dof = 2*sum(n_SJ3>0) - npars;
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
if Plot, PlotIt(pars, LevSJ3, ChosenModel, SJ3_TF, SJ3_S, SJ3_RF); end
% rearrange parameters to return
params = [pars(1:4) epskap(1:9)];
if ChosenModel==0, params(5:13)=NaN;
elseif ChosenModel==2, params(7)=NaN; params(12:13)=NaN;
elseif ChosenModel==3, params(6)=NaN; params(10:11)=NaN;
elseif ChosenModel==4, params(5)=NaN; params(8:9)=NaN;
elseif ChosenModel==5, params(6:7)=NaN; params(10:13)=NaN;
elseif ChosenModel==6, params(5)=NaN; params(7)=NaN; params(8:9)=NaN; params(12:13)=NaN;
elseif ChosenModel==7, params(5:6)=NaN; params(8:11)=NaN;
end
% create report
report = struct('Problem', ['Separate fit of SJ3 data with ' output.algorithm],...
                'FlagFrom_fmincon', Final_flag, ...
                'NumIterations', Final_iterations, ...
                'NumFuncCount', Final_funcCount, ...
                'UserModel', Model, ...
                'FittedModel', ChosenModel, ...
                'NumFreeParameters', npars, ...
                'NumCells', 3*sum(n_SJ3>0), ...
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
                'SJ3_LamT', params(1), ...
                'SJ3_LamR', params(2), ...
                'SJ3_Tau', params(3), ...
                'SJ3_Delta', params(4), ...
                'SJ3_epsilonTF', params(5), ...
                'SJ3_epsilonS', params(6), ...
                'SJ3_epsilonRF', params(7), ...
                'SJ3_kappa_TFintoS', params(8), ...
                'SJ3_kappa_TFintoRF', params(9), ...
                'SJ3_kappa_SintoTF', params(10), ...
                'SJ3_kappa_SintoRF', params(11), ...
                'SJ3_kappa_RFintoTF', params(12), ...
                'SJ3_kappa_RFintoS', params(13));
warning(s1); warning(s2);
% -----------------Calculate Likelihood------------------------------------
function L = LogLikelihood(x, SJ3_TF, SJ3_S, SJ3_RF, LevSJ3, Model)
  [lam_T lam_R tau delta_SJ3 epskap] = GetParVec(x, Model);
  penalty = 1.0e20;
  if lam_T<=0 || lam_R<=0 || delta_SJ3<0 || any(epskap<0) || any(epskap>1)
    L = penalty;
  else
    [P_SJ3_TF P_SJ3_S] = psi_SJ3_TF_S(LevSJ3, lam_T, lam_R, tau, delta_SJ3, epskap);
    P_SJ3_RF = 1 - P_SJ3_TF - P_SJ3_S;
    L = -nansum([SJ3_TF.*log(P_SJ3_TF) SJ3_S.*log(P_SJ3_S) SJ3_RF.*log(P_SJ3_RF)]);
    if ~isreal(L) || isinf(L), L = penalty; end
  end
% -----------------Psychometric Functions
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
function [lam_T lam_R tau delta_SJ3 epskap] = GetParVec(p, Model)
lam_T = p(1);  lam_R = p(2);  tau = p(3); delta_SJ3 = p(4);
if Model==0, epskap = zeros(1,9);
elseif Model==1, epskap = [p(5:8) 1-p(8) p(9) 1-p(9) p(10) 1-p(10)];
elseif Model==2, epskap = [p(5:6) 0 p(7) 1-p(7) p(8) 1-p(8) 0 0];
elseif Model==3, epskap = [p(5) 0 p(6) p(7) 1-p(7) 0 0 p(8) 1-p(8)];
elseif Model==4, epskap = [0 p(5) p(6) 0 0 p(7) 1-p(7) p(8) 1-p(8)];
elseif Model==5, epskap = [p(5) 0 0 p(6) 1-p(6) 0 0 0 0];
elseif Model==6, epskap = [0 p(5) 0 0 0 p(6) 1-p(6) 0 0];
elseif Model==7, epskap = [0 0 p(5) 0 0 0 0 p(6) 1-p(6)];
end
%------------------Distribution of the arrival-time difference-------------
function d = difcdf(x,shift,lam1,lam2)
  y = x-shift;
  left = min(lam1.*exp(lam2.*y)./(lam1+lam2),1);
  right = max(1-(lam2.*exp(-lam1.*y))./(lam1+lam2),0);
  d = (y<=0).*left + (y>0).*right;
% -----------------Plot Results--------------------------------------------
function PlotIt(pars, LevSJ3, Model, SJ3_TF, SJ3_S, SJ3_RF)
[lam_T lam_R tau delta_SJ3 epskap] = GetParVec(pars, Model);
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/10 scrsz(4)/3 4*scrsz(3)/15 scrsz(4)/4])
lw = 3; % line width
ms = 8; % marker size
fsa = 12; % font size for axis labels and title
FS = ['\fontsize{',num2str(fsa),'}'];
fsl = 11; % font size for legend
  subplot (1,1,1)
xinf = min(LevSJ3);  xsup = max(LevSJ3);
step = (xsup-xinf)/800;
xval = xinf:step:xsup; % values x-axis
[P_SJ3_TF P_SJ3_S] = psi_SJ3_TF_S(xval, lam_T, lam_R, tau, delta_SJ3, epskap);
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
