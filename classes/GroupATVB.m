classdef GroupATVB < Group
      
    properties
        latencies = 13
    end
    
    methods
       
        function self = GroupATVB(len, xlabels, xdata, titleLabels, ylimits)
           self@Group(len, xlabels, xdata, titleLabels, ylimits);
        end
               
        function plotSubjectsGFit(self, varargin)
            
            [filtered_subjects, subjs_title] = self.filterSubjects(varargin{:});

            [stat_a_tv, stat_t_av, stat_v_at] = self.getSubjectsStat(filtered_subjects);
            
            stat_a_tv.gfit = StatsFits.gaussianFit(stat_a_tv.data_mean, self.xdata);
            stat_t_av.gfit = StatsFits.gaussianFit(stat_t_av.data_mean, self.xdata);
            stat_v_at.gfit = StatsFits.gaussianFit(stat_v_at.data_mean, self.xdata);
            
            title_a_tv = [subjs_title " A vs TV"];
            title_t_av = [subjs_title " T vs AV"];
            title_v_at = [subjs_title " V vs AT"];            
            
            self.plotDataWithError(stat_a_tv.data_mean, stat_a_tv.data_sem, stat_a_tv.gfit.mu, stat_a_tv.gfit.sigma, title_a_tv);
            self.plotDataWithError(stat_t_av.data_mean, stat_t_av.data_sem, stat_t_av.gfit.mu, stat_t_av.gfit.sigma, title_t_av);
            self.plotDataWithError(stat_v_at.data_mean, stat_v_at.data_sem, stat_v_at.gfit.mu, stat_v_at.gfit.sigma, title_v_at);
            
%             self.plotDataWithErrorFit(stat_a_tv.data_mean, stat_a_tv.data_sem, stat_a_tv.fit.y, stat_a_tv.fit.mu, stat_a_tv.fit.sigma, title_a_tv);
%             self.plotDataWithErrorFit(stat_t_av.data_mean, stat_t_av.data_sem, stat_t_av.fit.y, stat_t_av.fit.mu, stat_t_av.fit.sigma, title_a_tv);
%             self.plotDataWithErrorFit(stat_v_at.data_mean, stat_v_at.data_sem, stat_v_at.fit.y, stat_v_at.fit.mu, stat_v_at.fit.sigma, title_a_tv);
        end
                
        function create_tabbed_data(self, filename)
            
            fid = fopen(filename, 'w');

            fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                         'label','group','age','gender','experiment','task','ntrials', 'err_sim', 'err_1_23', 'err_23_1', 'mu','sigma', ...
                         'TFsb', 'RFsb', 'sr', 'mp', 'peak');

            for k=1:self.number
                self.subjects{k}.writeData(fid);
            end
            fclose(fid);

        end 
        
        % returns 3 struct(data_mean, data_sd, data_sem)
        function [stat_a_tv, stat_t_av, stat_v_at] = getSubjectsStat(self, filtered_subjects)
            
            % CALCULATE DATA
            len = length(filtered_subjects);
            
            data_a_tv = zeros(len, self.latencies);
            data_t_av = zeros(len, self.latencies);
            data_v_at = zeros(len, self.latencies);
            
            for s=1:length(filtered_subjects)
                
                subj = filtered_subjects{s};
                data_a_tv(s,:) = subj.a_tv.data;
                data_t_av(s,:) = subj.t_av.data;
                data_v_at(s,:) = subj.v_at.data;
            end
            
            stat_a_tv = StatsFits.getStat(data_a_tv);
            stat_t_av = StatsFits.getStat(data_t_av);
            stat_v_at = StatsFits.getStat(data_v_at);
           
        end   
        
        function [stat_a_tv, stat_t_av, stat_v_at] = plotSubjectsSJ2(self, varargin)
            
            [filtered_subjects, ~]              = self.filterSubjects(varargin{:});
            [stat_a_tv, stat_t_av, stat_v_at]   = self.getSubjectsStat(filtered_subjects);  % returns 3 struct(data_mean, data_sd, data_sem)
            
            data_a_tv = zeros(3, self.latencies);
            data_t_av = zeros(3, self.latencies);
            data_v_at = zeros(3, self.latencies);
            
            data_a_tv(1,:) = self.xdata;
            data_a_tv(2,:) = (100 - stat_a_tv.data_mean);            
            data_a_tv(3,:) = stat_a_tv.data_mean;            
            
            data_t_av(1,:) = self.xdata;
            data_t_av(2,:) = (100 - stat_t_av.data_mean);            
            data_t_av(3,:) = stat_t_av.data_mean;            
            
            data_v_at(1,:) = self.xdata;
            data_v_at(2,:) = (100 - stat_v_at.data_mean);            
            data_v_at(3,:) = stat_v_at.data_mean;            
            
            stat_a_tv = StatsFits.sj2Fit(data_a_tv, 'A vs TV');
            stat_t_av = StatsFits.sj2Fit(data_t_av, 'T vs AV');
            stat_v_at = StatsFits.sj2Fit(data_v_at, 'V vs AT');
        end
   end
end