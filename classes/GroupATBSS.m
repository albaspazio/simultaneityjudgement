classdef GroupATBSS < Group
      
    properties
        latencies = 15
    end
    
    methods
       
        function self = GroupATBSS(len, xlabels, xdata, titleLabels, ylimits)
           self@Group(len, xlabels, xdata, titleLabels, ylimits);
        end
    
        
        function stat_a_t = getSubjectsStat(self, filtered_subjects)
            
            % CALCULATE DATA
            len = length(filtered_subjects);
            
            data_a_t = zeros(len, self.latencies);
            
            for s=1:length(filtered_subjects)
                
                subj = filtered_subjects{s};
                data_a_t(s,:) = subj.a_t.data;
            end
            
            stat_a_t = StatsFits.getStat(data_a_t, self.xdata);
        end        
        
        
        function plotSubjects(self, varargin)

            [filtered_subjects, subjs_title] = self.filterSubjects(varargin{:});

            stat_a_t  = self.getSubjectsStat(filtered_subjects);                        
            title_a_t = [subjs_title " A vs TV"];
            
            self.plotDataWithError(stat_a_t.data_mean, stat_a_t.data_sem, stat_a_t.fit.mu, stat_a_t.fit.sigma, title_a_t);
            
            
%             data_a_t = zeros(1,self.latencies);
            
%             if(isempty(varargin))
%                 % average all subjects
%                 filtered_subjects = self.subjects{:};
%                 subjs_title = "All subjects:";
%             else
%                 
%             end
%             
%             for s=1:self.number
%                 subj = self.subjects{s};
%                 data_a_t   = data_a_t + subj.a_t.data;
%             end
%             
%             data_a_t    = data_a_t./self.number;
%             fit_a_t     = StatsFits.gaussianFit(data_a_t, self.xdata);   % calculate group fit
%             title_a_t   = [subjs_title " A vs T"];
%             self.plotCurve(fit_a_t.y, data_a_t, fit_a_t.mu, fit_a_t.sigma, title_a_t);
        end
   end
end