classdef StatsFits

    
    properties
        
    end
    
    methods (Static)
        
        % data is [nsubj, ndelays] 
        function stat = getStat(data)  
            
            subjlen         = size(data,1);
            
            stat.data_mean  = mean(data, 1);
            stat.data_std   = std(data, 0, 1);
            stat.data_sem   = stat.data_std/sqrt(subjlen);
        end
        
        function gf = gaussianFit(data, xdata)

            try
                datafit = fit(xdata.', data.','gauss1');
                
                gf.sigma = datafit.c1/2;
                gf.mu    = datafit.b1;
                gf.y     = datafit.a1*exp(-((xdata-datafit.b1).^2)/datafit.c1^2);
            catch err
                gf.sigma = NaN;
                gf.mu    = NaN;
                gf.y     = zeros(length(data),1);        
            end
        end     
        
        % data[3,lat]
        function datafit = sj2Fit(data, title)
            
            LamBounds   = [1/200 1/3]; 
            TauBounds   = [-Inf Inf]; 
            DeltaBounds = [0 Inf];
            LamTStart   = [1/70 1/10]; 
            LamRStart   = [1/70 1/10]; 
            TauStart    = [-70 70];
            DeltaStart  = [20 150]; 
            ErrStart    = [0.05]; 
            Model       = 1; 
            Disp        = true;        

            if isempty(title)
                Plot        = false; 
                PlotRange   = [];
            else
                Plot        = true; 
                PlotRange   = [data(1, 1) data(1, end)];
            end

            datafit     = fit_SJ2(  data, LamBounds, TauBounds, DeltaBounds, ...
                                    LamTStart, LamRStart, TauStart, ...
                                    DeltaStart, ErrStart, Model, Plot, Disp, ['fit ' title]);
                                
%             SampleSize  = 1000;
%             ConfCoef    = 95;
%             FixedSeed   = true;
%             bstrap_data = bootstrap_SJ2(datafit, data, SampleSize, ConfCoef, FixedSeed);

            task = 'SJ2'; 
            Xi = 0;
            datafit.pm = PerformanceMeasures(task, datafit.SJ2_LamT, datafit.SJ2_LamR, datafit.SJ2_Tau, datafit.SJ2_Delta, Xi, PlotRange, ['perf ' title]);                                
            
            if strcmp(datafit.pm.SJ2_TFsb, 'Does not exist')
                datafit.pm.SJ2_TFsb = 0;
            end
            if strcmp(datafit.pm.SJ2_RFsb, 'Does not exist')
                datafit.pm.SJ2_RFsb = 0;
            end
            if strcmp(datafit.pm.SJ2_sr, 'Does not exist')
                datafit.pm.SJ2_sr = 0;
            end
            if strcmp(datafit.pm.SJ2_mp, 'Does not exist')
                datafit.pm.SJ2_mp = 0;
            end
            if strcmp(datafit.pm.SJ2_peak, 'Does not exist')
                datafit.pm.SJ2_peak = 0;
            end
            
            datafit.pm.SJ2_sr = datafit.pm.SJ2_sr/2;
        end
        
    end
end

