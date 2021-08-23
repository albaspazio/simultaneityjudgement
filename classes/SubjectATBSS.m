% it can manage either 13 or 15 latencies

classdef SubjectATBSS < Subject
    
    properties
        a_t
        errors
        err_control     % highest delay, control condition

    end
    
    methods
        function self = SubjectATBSS(varargin)
            self@Subject(varargin{1:6}); 
            
            if(nargin == 7)  % it includes also data to parse
                self = self.processData(varargin{7});
                self = self.getErrors();
                self = self.gfitData();
                self = self.SJ2fitData();
            end
        end
            
        function self = processData(self, data)
            
            data_a_t        = zeros(1, self.latencies);
            cnt_data_a_t    = zeros(1, self.latencies);

            midlatency      = ceil(self.latencies/2);  ...13 --> 7, 15 --> 8

            for n=1:length(data.type)

                ... determine value to insert
                if(strcmpi(data.answer(n,:), 'si') == true)
                    value = 1;
                else
                    value = 0;
                end

                if(data.type(n) == 12)

                    ... put simultaneous trials in all three arrays
                    data_a_t(midlatency)        = data_a_t(midlatency) + value;
                    cnt_data_a_t(midlatency)    = cnt_data_a_t(midlatency) + 1;

                else
                    ... shifted trials
                    switch(num2str(data.type(n)))
                        case '102'
                            first_smaller = true;
                            [data_a_t, cnt_data_a_t] = self.putValue(value, data_a_t, cnt_data_a_t, data.delay(n), first_smaller);

                        case '201'
                            first_smaller = false;
                            [data_a_t, cnt_data_a_t] = self.putValue(value, data_a_t, cnt_data_a_t, data.delay(n), first_smaller);
                    end
                end
            end
            
            self.a_t.data    = (data_a_t*100)./cnt_data_a_t;
            self.a_t.ntrial  = sum(cnt_data_a_t);            
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
            self.a_t.gfit = StatsFits.gaussianFit(self.a_t.data, self.xdata);
        end      
                
        function self = SJ2fitData(self)
            
            self.a_t.sj2fit = StatsFits.sj2Fit(self.getSJ2Data(), '');
        end 
        
        % returns [3,3,lat]
        function data = getSJ2Data(self)
            
            data = zeros(3, self.latencies);
            data(1,:) = self.xdata;
            data(2,:) = (100 - self.a_t.data);
            data(3,:) = self.a_t.data;
        end
        
        function self = getErrors(self)
            
            midlatency = ceil(self.latencies/2);  ...13 --> 7, 15 --> 8
            
            for t=1:length(self.a_t.data)
                if t == midlatency
                    self.errors(t) = 100 - self.a_t.data(t);
                else
                    self.errors(t) = self.a_t.data(t);
                end
            end
            
            self.err_sim = self.errors(midlatency);
            self.err_control = mean([self.errors(1) self.errors(end)]);
            
        end
        
        function writeData(self, fid)
            
            fprintf(fid, '%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n', ...
                         self.label, self.group, self.age, self.gender, self.experiment, 'A_T', self.a_t.ntrial, ...
                         self.a_t.gfit.mu, self.a_t.gfit.sigma,self.a_t.sj2fit.pm.SJ2_TFsb, self.a_t.sj2fit.pm.SJ2_RFsb, self.a_t.sj2fit.pm.SJ2_sr, self.a_t.sj2fit.pm.SJ2_mp, self.a_t.sj2fit.pm.SJ2_peak, ...
                         self.err_sim, self.err_control);
        end
        
       
        function writeErrorsData(self, fid)
            
            for l=1:ceil(self.latencies/2)
                fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%4.2f\n', ...
                             self.label, self.group, self.age, self.gender, self.experiment, 'A_T', self.a_t.ntrial, l, mean([self.errors(l) self.errors(self.latencies + 1 - l)]));
            
%             for l=1:self.latencies
%                 fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%4.2f\n', ...
%                              self.label, self.group, self.age, self.gender, self.experiment, 'A_T', self.a_t.ntrial, l, self.errors(l));
            end
        end
        
        function plotData(self, xdata, titles, ylimits)
            
            d1 = self.a_t.data;
            y1 = self.a_t.gfit.y;

            figure
            hold on; 
            plot(xdata,y1,'*b');
            plot(xdata,d1,'.k');
            title([self.label titles{1}]);
            ylim(ylimits)
            text(6,20, ['mu = ' num2str(self.a_t.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.a_t.gfit.sigma)]);
                
        end     ... end plotData
            
    end    ... end methods
end

