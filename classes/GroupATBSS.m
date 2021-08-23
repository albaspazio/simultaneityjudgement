classdef GroupATBSS < Group
    
    methods
       
        function self = GroupATBSS(varargin) % len, xlabels, xdata, titleLabels, ylimits
           self@Group(varargin{:});
        end
    
        function stat_a_t = getSubjectsStat(self, filtered_subjects)
            
            % CALCULATE DATA
            len = length(filtered_subjects);
            
            data_a_t = zeros(len, self.latencies);
            
            for s=1:length(filtered_subjects)
                
                subj = filtered_subjects{s};
                data_a_t(s,:) = subj.a_t.data;
            end
            
            stat_a_t = StatsFits.getStat(data_a_t);
        end        
        
        
            
        function create_tabbed_data(self, filename, errors_file)
            
            %-----------------------------------------------------
            % fits
            fid = fopen(filename, 'w');
            fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                         'subj','group','age','gender','experiment','task','ntrials', 'mu', 'sigma', 'TFsb', 'RFsb', 'sr', 'mp', 'peak', 'err_sim', 'err_control');

            for k=1:self.number
                self.subjects{k}.writeData(fid);
            end
            fclose(fid);
            
            %-----------------------------------------------------
            % errors
            feid = fopen(errors_file, 'w');
            fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                         'subj','group','age','gender','experiment','task','ntrials', 'latency', 'error');

            for k=1:self.number
                self.subjects{k}.writeErrorsData(feid);
            end
            
            fclose(feid);
        end          
        
        %% PLOT DATA  
        function plotSubjectsGFit(self, varargin)

            [filtered_subjects, subjs_title]    = self.filterSubjects(varargin{:});
            stat_a_t                            = self.getSubjectsStat(filtered_subjects);
            stat_a_t.gfit                       = StatsFits.gaussianFit(stat_a_t.data_mean, self.xdata);
            title_a_t                           = [subjs_title " A vs T"];
            
            self.plotDataWithError(stat_a_t.data_mean, stat_a_t.data_sem, stat_a_t.gfit.mu, stat_a_t.gfit.sigma, title_a_t);
            
        end
      
        function stat_a_t = plotSubjectsSJ2(self, varargin)
            
            [filtered_subjects, ~]  = self.filterSubjects(varargin{:});
            stat_a_t                = self.getSubjectsStat(filtered_subjects);  % returns one struct(data_mean, data_sd, data_sem)
            
            data_a_t = zeros(3, self.latencies);
            
            data_a_t(1,:) = self.xdata;
            data_a_t(2,:) = (100 - stat_a_t.data_mean);            
            data_a_t(3,:) = stat_a_t.data_mean;            
            
            
            stat_a_t = StatsFits.sj2Fit(data_a_t, 'A vs T');
        end        
   end
end