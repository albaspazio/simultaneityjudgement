SJ2data     = [-300 -240 -180 -120 -60 0 60 120 180 240 300; ...
                127 117 109 83 27 18 10 29 74 85 108; ...
                  3 3 11 27 83 112 90 101 36 15 12];
LamBounds   = [1/200 1/3]; 
TauBounds   = [-Inf Inf]; 
DeltaBounds = [0 Inf];
LamTStart   = [1/70 1/10]; 
LamRStart   = [1/70 1/10]; 
TauStart    = [-70 70];
DeltaStart  = [20 150]; 
ErrStart    = [0.05]; 
Model = 1; Plot = true; Disp = true;
%BiasStart = [0.5];  ...unused in SJ2



o_SJ2 = fit_SJ2(SJ2data, LamBounds, TauBounds, DeltaBounds, ...
                LamTStart, LamRStart, TauStart, ...
                DeltaStart, ErrStart, Model, Plot, Disp);

            
            
SampleSize = 100;
ConfCoef = 95;
FixedSeed = true;
bstrap_SJ2 = bootstrap_SJ2(o_SJ2, SJ2data, SampleSize, ConfCoef, FixedSeed);


task = 'SJ2'; 
PlotRange = [-400 400];
LamT = 0.026334; 
LamR = 0.021859; 
Tau = -28.062; 
Delta = 126.63; 
Xi = 0;
SJ2pm = PerformanceMeasures(task, LamT, LamR, Tau, Delta, Xi, PlotRange);


