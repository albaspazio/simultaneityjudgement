 classdef SubjectATVB < Subject
    
    properties
        a_tv
        t_av
        v_at
        latencies = 13
        err_simult
    end
    
    methods
        
        function self = SubjectATVB(varargin)
            self@Subject(varargin{1:6}); 
            
            if(nargin == 7)  % it includes also data to parse
                self = self.processData(varargin{7});
                self = self.getErrors();
                self = self.gfitData();
                self = self.SJ2fitData();
            end
        end            
            
       
        function self = processData(self, data)
            
            data_a_tv = zeros(1,self.latencies);
            data_t_av = zeros(1,self.latencies);
            data_v_at = zeros(1,self.latencies);

            cnt_data_a_tv = zeros(1,self.latencies);
            cnt_data_t_av = zeros(1,self.latencies);
            cnt_data_v_at = zeros(1,self.latencies);

            for n=1:length(data.type)

                ... determine value to insert
                if(strcmpi(data.answer(n,:), 'si') == true)
                    value = 1;
                else
                    value = 0;
                end

                if(data.type(n) == 0)

                    ... put simultaneous trials in all three arrays
                    data_a_tv(7) = data_a_tv(7) + value;
                    data_t_av(7) = data_t_av(7) + value;
                    data_v_at(7) = data_v_at(7) + value;

                    cnt_data_a_tv(7) = cnt_data_a_tv(7) + 1;
                    cnt_data_t_av(7) = cnt_data_t_av(7) + 1;
                    cnt_data_v_at(7) = cnt_data_v_at(7) + 1;

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
            self.a_tv.ntrial = sum(cnt_data_a_tv);
            self.t_av.data   = data_t_av;
            self.t_av.ntrial = sum(cnt_data_t_av);
            self.v_at.data   = data_v_at;
            self.v_at.ntrial = sum(cnt_data_v_at);
            
            self.ntrials     = length(data.type);
            
            self.err_simult  = 100 - data_a_tv(7);
        end
        
        function [arr, cnt_arr] = putValue(self, val, arr, cnt_arr, delay, first_smaller)

            id = 1;
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
            
            self.err_simult     = 100 - self.a_tv.data(7);
            
            self.a_tv.err_a_tv  = 100 - mean(self.a_tv.data(1:6));
            self.a_tv.err_tv_a  = 100 - mean(self.a_tv.data(8:13));
            
            self.t_av.err_t_av  = 100 - mean(self.t_av.data(1:6));
            self.t_av.err_av_t  = 100 - mean(self.t_av.data(8:13));
            
            self.v_at.err_v_at  = 100 - mean(self.v_at.data(1:6));
            self.v_at.err_at_v  = 100 - mean(self.v_at.data(8:13));            

        end
        
        % label, group, age, gender, experiment, task, ntrials,
        % err_simult, err_1_23, err_23_1, mu, sigma,
        % TFsb, RFsb, sr, mp, peak
        function writeData(self, fid)
            
            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'A_TV', self.a_tv.ntrial, ...
                         self.err_simult, self.a_tv.err_a_tv, self.a_tv.err_tv_a, self.a_tv.gfit.mu, self.a_tv.gfit.sigma,...
                         self.a_tv.sj2fit.pm.SJ2_TFsb, self.a_tv.sj2fit.pm.SJ2_RFsb, self.a_tv.sj2fit.pm.SJ2_sr, self.a_tv.sj2fit.pm.SJ2_mp, self.a_tv.sj2fit.pm.SJ2_peak);

            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'T_AV', self.t_av.ntrial, ...
                         self.err_simult, self.t_av.err_t_av, self.t_av.err_av_t, self.t_av.gfit.mu, self.t_av.gfit.sigma,...
                         self.t_av.sj2fit.pm.SJ2_TFsb, self.t_av.sj2fit.pm.SJ2_RFsb, self.t_av.sj2fit.pm.SJ2_sr, self.t_av.sj2fit.pm.SJ2_mp, self.t_av.sj2fit.pm.SJ2_peak);

            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'V_AT', self.v_at.ntrial, ...
                         self.err_simult, self.v_at.err_v_at, self.v_at.err_at_v, self.v_at.gfit.mu, self.v_at.gfit.sigma,...
                         self.v_at.sj2fit.pm.SJ2_TFsb, self.v_at.sj2fit.pm.SJ2_RFsb, self.v_at.sj2fit.pm.SJ2_sr, self.v_at.sj2fit.pm.SJ2_mp, self.v_at.sj2fit.pm.SJ2_peak);
        end
        
        function plotData(self, xdata, titles)
            
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
            text(6,20, ['mu = ' num2str(self.a_tv.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.a_tv.gfit.sigma)]);

            figure
            hold on; 
            plot(xdata,y2,'*b');
            plot(xdata,d2,'.k');
            title([self.label titles{2}]);
            text(6,20, ['mu = ' num2str(self.t_av.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.t_av.gfit.sigma)]);

            figure
            hold on; 
            plot(xdata,y3,'*b');
            plot(xdata,d3,'.k');
            title([self.label titles{3}]);
            text(6,20, ['mu = ' num2str(self.v_at.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.v_at.gfit.sigma)]);
        end     ... end plotData
            
    end    ... end methods
end

