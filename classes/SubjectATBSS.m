classdef SubjectATBSS < Subject
    
    properties
        a_t
    end
    
    methods
        function self = SubjectATBSS(lab, ag, gend, grp, exp)
            self@Subject(lab, ag, gend, grp, exp); 
        end            
            
        function self = processData(self, data)
            
            data_a_t = zeros(1,15);
            cnt_data_a_t = zeros(1,15);
    
            for n=1:length(data.type)

                ... determine value to insert
                if(strcmpi(data.answer(n,:), 'si') == true)
                    value = 1;
                else
                    value = 0;
                end

                if(data.type(n) == 0)

                    ... put simultaneous trials in all three arrays
                    data_a_t(8)        = data_a_t(8) + value;
                    cnt_data_a_t(8)    = cnt_data_a_t(8) + 1;

                else
                    ... shifted trials
                    switch(num2str(data.type(n)))
                        case '3'
                            first_smaller = true;
                            [data_a_t, cnt_data_a_t] = self.putValue(value, data_a_t, cnt_data_a_t, data.delay(n), first_smaller);

                        case '4'
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
            arr(id) = arr(id) + val;
            cnt_arr(id) = cnt_arr(id) + 1;    
        end        
        
        function self = setFit(self, xdata)

            self.a_t.gfit = StatsFits.gaussianFit(self.a_t.data, xdata);
 
        end      
        
        function writeData(self, fid)
            
            fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%4.2f\t%4.2f\n', self.label, self.group, self.age, self.gender, self.experiment, ...
                                                                     'A_T', self.a_t.ntrial, self.a_t.gfit.mu, self.a_t.gfit.sigma);
        end
        
        function plotData(self, xdata, titles)
            
            d1 = self.a_t.data;
            y1 = self.a_t.fit.y;

            figure
            hold on; 
            plot(xdata,y1,'*b');
            plot(xdata,d1,'.k');
            title([self.label titles{1}]);
            text(6,20, ['mu = ' num2str(self.a_t.gfit.mu)]);
            text(6,10, ['sigma = ' num2str(self.a_t.gfit.sigma)]);
                
        end     ... end plotData
            
    end    ... end methods
end

