classdef Group
    properties
        subjects;
        number = 0
        xlabels
        xdata
        titleLabels
        ylimits
        latencies

    end
      
    methods
       
        function self = Group(len, xlabels, xdata, titleLabels, ylimits)
           self.subjects    = cell(len,1);
           self.number      = 0;
           self.xlabels     = xlabels;
           self.xdata       = xdata;
           self.titleLabels = titleLabels;
           self.ylimits     = ylimits;
           self.latencies   = length(xdata);
        end
    
        function self = add(self, subj)
              self.number                 = self.number + 1;
              self.subjects{self.number}  = subj;
        end
        
        %% GET/FILTER SUBJECTS
        % returns subset of subjects' cellarray 
        function [filtered_subjects, subjs_title] = filterSubjects(self, varargin)
            
            if(isempty(varargin))
                % average all subjects
                filtered_subjects = self.subjects(:);
                subjs_title = "All subjects:";
            else
                % pairs of : field | values = [] are present. 
                % e.g. "age", [7,8]
                subjs_title = Group.filterString(varargin{:});
                
                % extract first pair
                field = varargin{1};
                values = varargin{2};
                indexes = self.extractSubset(field, values);                
                
                % then I intersect it with following ones
                nvars = length(varargin);
                for v=2:nvars/2
                   field = varargin{2*v-1};
                   values = varargin{2*v};
                   indexes = intersect(indexes, self.extractSubset(field, values));
                end
                filtered_subjects = self.subjects(indexes);
            end            
        end
        
        % assumes values is an array with values of the same type
        function indexes = extractSubset(self, field, values)

            fieldvalues = cell(1, self.number);
            isnum       = isnumeric(self.subjects{1}.(field));
                    
            for s=1:self.number
                if isnum == true
                    fieldvalues{s} = num2str(self.subjects{s}.(field));
                else
                    fieldvalues{s} = self.subjects{s}.(field);
                end
            end
                
            indexes = [];
            for v=1:length(values)
                v_ids = find(strcmp(fieldvalues, num2str(values(v))) == true);
                indexes = [indexes v_ids];
            end
        end

        % returns a subject instance
        function subj = getSubjectByLabel(self, label)
            
            for s=1:self.number
               if strcmp(self.subjects{s}.label, label)
                   subj = self.subjects{s};
                   break;
               end
            end
        end
        
        %% PLOT DATA
        function plotSubject(self, label, xdata, titleLabels)
            self.getSubjectByLabel(label).plotData(xdata, titleLabels, self.ylimits);
        end
        
        function plotDataWithError(self, data, datastd, mu, sigma, titl)
   
            figure
            hold on; 
            
            errorbar(self.xdata, data, datastd, 'Color', 'b');
            title(titl);
            plot(self.xdata, data, 'Marker', 'o', 'MarkerFaceColor', 'b');
            ylim(self.ylimits);
            text(0,20, ['mu = ' num2str(mu)]);
            text(0,17, ['sigma = ' num2str(sigma)]);         
        end
        
        function plotDataWithErrorFit(self, data, datastd, fity, mu, sigma, titl)
   
            figure
            hold on; 
            
            errorbar(self.xdata, data, datastd);
            title(titl);
            plot(self.xdata, fity, 'Marker', '*'); ..., 'MarkerFaceColor', 'b'); ..., 'LineStyle', '-', 'Color', 'b');
            ylim(self.ylimits);
            text(6,20, ['mu = ' num2str(mu)]);
            text(6,10, ['sigma = ' num2str(sigma)]);         
        end
        
        function plotCurve(self, y, d, mu, sigma, titl)

            figure
            hold on; 
            ...plot(self.xdata, y, '*b', '-x');
            plot(self.xdata, y, 'Marker', '*', 'MarkerFaceColor', 'b', 'LineStyle', '-', 'Color', 'b');
            plot(self.xdata, d, 'Marker', 'x', 'MarkerFaceColor', 'g', 'LineStyle', '--', 'Color', 'g');
            title(titl);
            ylim(self.ylimits)
            text(6,20, ['mu = ' num2str(mu)]);
            text(6,10, ['sigma = ' num2str(sigma)]);
        end
    end
      
   methods(Static)
       
        function str = filterString(varargin)
            nvars = length(varargin);
            str = "";
            for v=1:nvars/2
               field = varargin{2*v-1};
               values = varargin{2*v};
               str = strcat(str, field, "=", mat2str(values), "; ");
            end
        end      
      
   end    
    
end