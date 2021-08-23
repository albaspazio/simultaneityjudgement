 classdef SubjectATVB < Subject
    
    properties
        a_tv
        t_av
        v_at
        errors_a_tv
        errors_t_av
        errors_v_at

        err_a_tv_control     % highest delay, control condition
        err_t_av_control     % highest delay, control condition
        err_v_at_control     % highest delay, control condition

    end
    
    methods
        
        %% constructor
        function self = SubjectATVB(varargin)
            self@Subject(varargin{1:6}); 
            
            if(nargin == 7)  % it includes also data to parse
                self = self.processData(varargin{7});
                self = self.getErrors();
                self = self.gfitData();
                self = self.SJ2fitData();
            end
        end
        
        %%
        function self = processData(self, data)
            
            data_a_tv = zeros(1, self.latencies);
            data_t_av = zeros(1, self.latencies);
            data_v_at = zeros(1, self.latencies);

            cnt_data_a_tv = zeros(1, self.latencies);
            cnt_data_t_av = zeros(1, self.latencies);
            cnt_data_v_at = zeros(1, self.latencies);

            midlatency = ceil(self.latencies/2);  ...13 --> 7, 15 --> 8
            
            for n=1:length(data.type)

                ... determine value to insert
                if(strcmpi(data.answer(n,:), 'si') == true)
                    value = 1;
                else
                    value = 0;
                end

                if(data.type(n) == 0)

                    ... put simultaneous trials in all three arrays
                    data_a_tv(midlatency) = data_a_tv(midlatency) + value;
                    data_t_av(midlatency) = data_t_av(midlatency) + value;
                    data_v_at(midlatency) = data_v_at(midlatency) + value;

                    cnt_data_a_tv(midlatency) = cnt_data_a_tv(midlatency) + 1;
                    cnt_data_t_av(midlatency) = cnt_data_t_av(midlatency) + 1;
                    cnt_data_v_at(midlatency) = cnt_data_v_at(midlatency) + 1;

                else
                    ... shifted trials
                    switch(num2str(data.type(n)))
                        case '1'
                            first_smaller = true;
                            [data_a_tv, cnt_data_a_tv] = self.putValue(value, data_a_tv, cnt_data_a_tv, data.delay(n), first_smaller);

                        case '2'
                            first_smaller = false;
                            [data_a_tv, cnt_data_a_tv] = self.putValue(value, data_a_tv, cnt_data_a_tv, data.delay(n), first_smaller);

                        case '3'
                            first_smaller = true;
                            [data_t_av, cnt_data_t_av] = self.putValue(value, data_t_av, cnt_data_t_av, data.delay(n), first_smaller);

                        case '4'
                            first_smaller = false;
                            [data_t_av, cnt_data_t_av] = self.putValue(value, data_t_av, cnt_data_t_av, data.delay(n), first_smaller);

                        case '5'
                            first_smaller = true;
                            [data_v_at, cnt_data_v_at] = self.putValue(value, data_v_at, cnt_data_v_at, data.delay(n), first_smaller);

                        case '6'
                            first_smaller = false;
                            [data_v_at, cnt_data_v_at] = self.putValue(value, data_v_at, cnt_data_v_at, data.delay(n), first_smaller);
                    end
                end
            end
            
            data_a_tv = (data_a_tv*100)./cnt_data_a_tv;
            data_t_av = (data_t_av*100)./cnt_data_t_av;
            data_v_at = (data_v_at*100)./cnt_data_v_at;

            self.a_tv.data   = data_a_tv;
            self.t_av.data   = data_t_av;
            self.v_at.data   = data_v_at;
            
            self.a_tv.ntrial = sum(cnt_data_a_tv);
            self.t_av.ntrial = sum(cnt_data_t_av);
            self.v_at.ntrial = sum(cnt_data_v_at);
            
            self.ntrials     = length(data.type);
        end
        
        function [arr, cnt_arr] = putValue(self, val, arr, cnt_arr, delay, first_smaller)

            id = 1;
            if self.latencies == 13

                switch(num2str(delay))
                    case '50'
                        if(first_smaller == true)
                            id = 6;
                        else
                            id = 8;
                        end
                    case '100'
                        if(first_smaller == true)
                            id = 5;
                        else
                            id = 9;
                        end            
                    case '200'
                        if(first_smaller == true)
                            id = 4;
                        else
                            id = 10;
                        end            
                    case '300'
                        if(first_smaller == true)
                            id = 3;
                        else
                            id = 11;
                        end            
                    case '400'
                        if(first_smaller == true)
                            id = 2;
                        else
                            id = 12;
                        end            
                    case '800'
                        if(first_smaller == true)
                            id = 1;
                        else
                            id = 13;
                        end            
                end
            else
                switch(num2str(delay))
                    case '50'
                        if(first_smaller == true)
                            id = 7;
                        else
                            id = 9;
                        end
                    case '100'
                        if(first_smaller == true)
                            id = 6;
                        else
                            id = 10;
                        end            
                    case '200'
                        if(first_smaller == true)
                            id = 5;
                        else
                            id = 11;
                        end            
                    case '300'
                        if(first_smaller == true)
                            id = 4;
                        else
                            id = 12;
                        end            
                    case '400'
                        if(first_smaller == true)
                            id = 3;
                        else
                            id = 13;
                        end            
                    case '800'
                        if(first_smaller == true)
                            id = 2;
                        else
                            id = 14;
                        end  
                    case '1200'
                        if(first_smaller == true)
                            id = 1;
                        else
                            id = 15;
                        end                          
                end                
            end
            
            arr(id)     = arr(id) + val;
            cnt_arr(id) = cnt_arr(id) + 1;     
        end        
        
        function self = gfitData(self)
            
            self.a_tv.gfit = StatsFits.gaussianFit(self.a_tv.data, self.xdata);
            self.t_av.gfit = StatsFits.gaussianFit(self.t_av.data, self.xdata);
            self.v_at.gfit = StatsFits.gaussianFit(self.v_at.data, self.xdata);
        end      
                
        function self = SJ2fitData(self)
            
            data = self.getSJ2Data();
            
            self.a_tv.sj2fit = StatsFits.sj2Fit(squeeze(data(1,:,:)), '');
            self.t_av.sj2fit = StatsFits.sj2Fit(squeeze(data(2,:,:)), '');
            self.v_at.sj2fit = StatsFits.sj2Fit(squeeze(data(3,:,:)), '');
        end 
        
        % returns [3,3,lat]
        function data = getSJ2Data(self)
            
            data = zeros(3, 3, self.latencies);

            data(1,1,:) = self.xdata;
            data(1,2,:) = (100 - self.a_tv.data);
            data(1,3,:) = self.a_tv.data;
            
            data(2,1,:) = self.xdata;
            data(2,2,:) = (100 - self.t_av.data);
            data(2,3,:) = self.t_av.data;
            
            data(3,1,:) = self.xdata;
            data(3,2,:) = (100 - self.v_at.data);
            data(3,3,:) = self.v_at.data;
        end
        
        function self = getErrors(self)
            
            midlatency = ceil(self.latencies/2);  ...13 --> 7, 15 --> 8
            
            for t=1:length(self.a_tv.data)
                if t == midlatency
                    self.errors_a_tv(t) = 100 - self.a_tv.data(t);
                else
                    self.errors_a_tv(t) = self.a_tv.data(t);
                end
            end
            
            for t=1:length(self.t_av.data)
                if t == midlatency
                    self.errors_t_av(t) = 100 - self.t_av.data(t);
                else
                    self.errors_t_av(t) = self.t_av.data(t);
                end
            end            
            
            for t=1:length(self.v_at.data)
                if t == midlatency
                    self.errors_v_at(t) = 100 - self.v_at.data(t);
                else
                    self.errors_v_at(t) = self.v_at.data(t);
                end
            end  
            
            self.err_sim    = self.errors_a_tv(midlatency);    % midlatency values are the same 
            
            self.err_a_tv_control = mean([self.errors_a_tv(1) self.errors_a_tv(end)]);            
            self.err_t_av_control = mean([self.errors_t_av(1) self.errors_t_av(end)]);            
            self.err_v_at_control = mean([self.errors_v_at(1) self.errors_v_at(end)]);            
        end
        
        % label, group, age, gender, experiment, task, ntrials,
        % mu, sigma, TFsb, RFsb, sr, mp, peak
        % err_sim, err_X_XX_control
        function writeData(self, fid)
            
            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'A_TV', self.a_tv.ntrial, ...
                         self.a_tv.gfit.mu, self.a_tv.gfit.sigma,...
                         self.a_tv.sj2fit.pm.SJ2_TFsb, self.a_tv.sj2fit.pm.SJ2_RFsb, self.a_tv.sj2fit.pm.SJ2_sr, self.a_tv.sj2fit.pm.SJ2_mp, self.a_tv.sj2fit.pm.SJ2_peak, ...
                         self.err_sim, self.err_a_tv_control);

            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'T_AV', self.t_av.ntrial, ...
                         self.t_av.gfit.mu, self.t_av.gfit.sigma,...
                         self.t_av.sj2fit.pm.SJ2_TFsb, self.t_av.sj2fit.pm.SJ2_RFsb, self.t_av.sj2fit.pm.SJ2_sr, self.t_av.sj2fit.pm.SJ2_mp, self.t_av.sj2fit.pm.SJ2_peak, ...
                         self.err_sim, self.err_t_av_control);

            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'V_AT', self.v_at.ntrial, ...
                         self.v_at.gfit.mu, self.v_at.gfit.sigma,...
                         self.v_at.sj2fit.pm.SJ2_TFsb, self.v_at.sj2fit.pm.SJ2_RFsb, self.v_at.sj2fit.pm.SJ2_sr, self.v_at.sj2fit.pm.SJ2_mp, self.v_at.sj2fit.pm.SJ2_peak, ...
                         self.err_sim, self.err_v_at_control);
        end
        
        % label, group, age, gender, experiment, task, ntrials, id latency, error
        function writeErrorsData(self, fid)
            
            for l=1:ceil(self.latencies/2)
                fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%4.2f\n', ...
                             self.label, self.group, self.age, self.gender, self.experiment, 'A_TV', self.a_tv.ntrial, l, mean([self.errors_a_tv(l) self.errors_a_tv(self.latencies + 1 - l)]));

                fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%4.2f\n', ...
                             self.label, self.group, self.age, self.gender, self.experiment, 'T_AV', self.t_av.ntrial, l, mean([self.errors_t_av(l) self.errors_t_av(self.latencies + 1 - l)]));

                fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%4.2f\n', ...
                             self.label, self.group, self.age, self.gender, self.experiment, 'V_AT', self.v_at.ntrial, l, mean([self.errors_v_at(l) self.errors_v_at(self.latencies + 1 - l)]));
            end
        end        
        
        function plotData(self, xdata, titles, ylimits)
            
            d1 = self.a_tv.data;
            d2 = self.t_av.data;
            d3 = self.v_at.data;

            y1 = self.a_tv.gfit.y;
            y2 = self.t_av.gfit.y;
            y3 = self.v_at.gfit.y;

            figure
            hold on; 
            plot(xdata,y1,'*b');
            plot(xdata,d1,'.k');
            title([self.label titles{1}]);
            ylim(ylimits)
            text(6,20, ['mu = ' num2str(self.a_tv.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.a_tv.gfit.sigma)]);

            figure
            hold on; 
            plot(xdata,y2,'*b');
            plot(xdata,d2,'.k');
            title([self.label titles{2}]);
            ylim(ylimits)
            text(6,20, ['mu = ' num2str(self.t_av.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.t_av.gfit.sigma)]);

            figure
            hold on; 
            plot(xdata,y3,'*b');
            plot(xdata,d3,'.k');
            title([self.label titles{3}]);
            ylim(ylimits)
            text(6,20, ['mu = ' num2str(self.v_at.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.v_at.gfit.sigma)]);
        end     ... end plotData
            
    end    ... end methods
end

