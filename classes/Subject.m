classdef Subject
    
    properties
        label
        age
        gender
        group
        experiment
        xdata
        ntrials
    end
    
    methods
        function obj = Subject(lab, ag, gend, grp, exp, xdata)
            obj.label       = lab;
            obj.age         = ag;
            obj.gender      = gend;
            obj.group       = grp;
            obj.experiment  = exp;
            obj.xdata       = xdata;
        end
    end
end

