classdef SubjectATVB < Subject
    
    properties
        a_tv
        t_av
        v_at
    end
    
    methods
        
        function self = SubjectATVB(lab, ag, gend, grp, exp)
            self@Subject(lab, ag, gend, grp, exp); 
        end            
            
        function self = processData(self, data)
            
            data_a_tv = zeros(1,13);
            data_t_av = zeros(1,13);
            data_v_at = zeros(1,13);

            cnt_data_a_tv = zeros(1,13);
            cnt_data_t_av = zeros(1,13);
            cnt_data_v_at = zeros(1,13);

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
        
        function self = setFit(self, xdata)

            try
                fit_a_tv = fit(xdata.',self.a_tv.data.','gauss1');
                self.a_tv.sigma = fit_a_tv.c1/2;
                self.a_tv.mu    = fit_a_tv.b1;
                self.a_tv.y     = fit_a_tv.a1*exp(-((xdata-fit_a_tv.b1).^2)/fit_a_tv.c1^2);
            catch err
                self.a_tv.sigma = NaN;
                self.a_tv.mu    = NaN;
                self.a_tv.y     = zeros(13,1);        
            end

            try
                fit_t_av = fit(xdata.',self.t_av.data.','gauss1');
                self.t_av.sigma = fit_t_av.c1/2;
                self.t_av.mu    = fit_t_av.b1;
                self.t_av.y     = fit_t_av.a1*exp(-((xdata-fit_t_av.b1).^2)/fit_t_av.c1^2);
            catch err
                self.t_av.sigma = NaN;
                self.t_av.mu    = NaN;
                self.t_av.y     = zeros(13,1);        
            end

            try
                fit_v_at = fit(xdata.',self.v_at.data.','gauss1');
                self.v_at.sigma = fit_v_at.c1/2;
                self.v_at.mu    = fit_v_at.b1;
                self.v_at.y     = fit_v_at.a1*exp(-((xdata-fit_v_at.b1).^2)/fit_v_at.c1^2);
            catch err
                self.v_at.sigma = NaN;
                self.v_at.mu    = NaN;
                self.v_at.y     = zeros(13,1);        
            end    
        end      
        
        function writeData(self, fid)
            
            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\n', self.label, self.group, self.age, self.gender, self.experiment, ...
                                                                     'A_TV', self.a_tv.ntrial, self.a_tv.mu, self.a_tv.sigma);

            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\n', self.label, self.group, self.age, self.gender, self.experiment, ...
                                                                     'T_AV', self.t_av.ntrial, self.t_av.mu, self.t_av.sigma);

            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\n', self.label, self.group, self.age, self.gender, self.experiment, ...
                                                                     'V_AT', self.v_at.ntrial, self.v_at.mu, self.v_at.sigma);             
        end
        
        function plotData(self, xdata, titles)
            
            d1 = self.a_tv.data;
            d2 = self.t_av.data;
            d3 = self.v_at.data;

            y1 = self.a_tv.y;
            y2 = self.t_av.y;
            y3 = self.v_at.y;

            figure
            hold on; 
            plot(xdata,y1,'*b');
            plot(xdata,d1,'.k');
            title([self.label titles{1}]);
            text(6,20, ['mu = ' num2str(self.a_tv.mu)]);
            text(6,10, ['sigma = ' num2str(self.a_tv.sigma)]);
            ...set(gca,'xtick',(1:13),'xticklabel') ... ,xlabels)            

            figure
            hold on; 
            plot(xdata,y2,'*b');
            plot(xdata,d2,'.k');
            title([self.label titles{2}]);
            ...set(gca,'xtick',(1:13),'xticklabel') ... ,xlabels)            
            text(6,20, ['mu = ' num2str(self.t_av.mu)]);
            text(6,10, ['sigma = ' num2str(self.t_av.sigma)]);

            figure
            hold on; 
            plot(xdata,y3,'*b');
            plot(xdata,d3,'.k');
            title([self.label titles{3}]);
            ...set(gca,'xtick',(1:13),'xticklabel') ... ,xlabels)            
            text(6,20, ['mu = ' num2str(self.v_at.mu)]);
            text(6,10, ['sigma = ' num2str(self.v_at.sigma)]);
        end     ... end plotData
            
    end    ... end methods
end

