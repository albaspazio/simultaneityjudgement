classdef Subject
    
    properties
        label
        age
        gender
        group
        experiment
        xdata
        ntrials
        latencies
        err_sim
    end
    
    methods
        function obj = Subject(lab, ag, gend, grp, exp, xdata)
            obj.label       = lab;
            obj.age         = ag;
            obj.gender      = gend;
            obj.group       = grp;
            obj.experiment  = exp;
            obj.xdata       = xdata;
            obj.latencies   = length(xdata);
        end
    end
end

