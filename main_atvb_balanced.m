
%% edit this
root_dir            = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/';
experiment_folder   = 'bambini_bolzaneto';   ... 'bambini_bolzaneto' 'adulti_iit_1'
task_folder         = 'ATVBSSB';
results_folder      = 'results';


%% init variables
xlabels         = {'-300', '-250', '-200','0', '200', '250','300'};
titleLabels     = {'T/V - A - V/T', 'A/T - V - T/A','A/V - T - V/A'};
ylimits         = [0, 100];
xdata           = [-300, -250, -200, 0, 200, 250, 300];


%% start processing
data_dir        = fullfile(root_dir, experiment_folder, task_folder);
result_file     = fullfile(pwd, results_folder, strcat(experiment_folder, "_", task_folder, ".dat"));

files           = dir(strcat(data_dir,'/*.txt'));
nsubj           = length(files);
subjects        = Group(nsubj);

%% calculate the three timeseries (A_TV, T_AV, V_AT) of mean-yes-responses for each of the 7 delays.
... e.g. [0, 25, 25, 75, 100, 87.5,  91.6,  87.5, 87.5, 87.5, 62.5, 50,25]
for f=1:nsubj
    file_name       = files(f).name; ...'agga_8_1_ATVBSSU_TD_18122020_112553.txt';
    data            = tdfread(fullfile(data_dir, file_name));
    filename_parts  = split(file_name, '_');

    subj            = SubjectATBSS(filename_parts{1}, str2double(filename_parts{2}), filename_parts{3}, filename_parts{5}, filename_parts{4});
    subj            = subj.processData(data);
    subj            = subj.setFit(xdata); % uses Curve Fitting toolbox
    
    subjects        = subjects.add(subj);
end

clear filename_parts file f files nsubj label

create_tabbed_data(result_file, subjects);

% =======================================================================================
% accessory functions
% =======================================================================================

function subjects = processSubject(subjects, s, label, data)

    data_atv        = zeros(3,7);
    cnt_data_atv    = zeros(3,7);
    
    for n=1:length(data.type)
        [data_atv, cnt_data_atv] = putValue(data.type(n), data.delay(n), data.answer(n,:), data_atv, cnt_data_atv);
    end
    
    
    data_atv = ((data_atv)*100)./cnt_data_atv;
    
    subjects(s).label = label;
    subjects(s).a_tv.data   = data_atv(1);
    subjects(s).a_tv.ntrial = sum(cnt_data_atv(1));
    subjects(s).t_av.data   = data_atv(2);
    subjects(s).t_av.ntrial = sum(cnt_data_atv(2));
    subjects(s).v_at.data   = data_atv(3);
    subjects(s).v_at.ntrial = sum(cnt_data_atv(3));
    
end

% type, delay, answer, data, cnt
function [arr, cnt_arr] = putValue(type, delay, answer, arr, cnt_arr)
    ...T_A_V  = 10, V_A_T  = 11, A_T_V  = 12, V_T_A  = 13, A_V_T  = 14, T_V_A  = 15

    % get id of the measure
    switch(num2str(delay))
        case '200'
            if(type == 10 || type == 12 || type == 14)
                id_delay = 3;
            else
                id_delay = 5;
            end  
            
        case '250'
            if(type == 10 || type == 12 || type == 14)
                id_delay = 2;
            else
                id_delay = 6;
            end            
            
        case '300'
            if(type == 10 || type == 12 || type == 14)
                id_delay = 1;
            else
                id_delay = 7;
            end            
        otherwise
            id_delay = 4;
            
    end    

    % set answer % id in data array
    value = 0;
    switch(type)
        case '10'
        case '11'
            if(strcmp(answer(:), 'audio'))
                value = 1;
            end 
            id_data = 1;
        case '12'
        case '13'
            if(strcmp(answer(:), 'tactile'))
                value = 1;
            end 
            id_data = 2;
        case '14'
        case '15'
            if(strcmp(answer(:), 'video'))
                value = 1;
            end             
            id_data = 3;
        otherwise
            value = 0;
            
    end
        
    arr(id_data, id_delay)      = arr(id_data, id_delay) + value;
    cnt_arr(id_data, id_delay)  = cnt_arr(id_data, id_delay) + 1;       

end

function subject_data = setFit(subject_data, xdata)
    
    try
        fit_a_tv = fit(xdata.',subject_data.a_tv.data.','gauss1');
        subject_data.a_tv.sigma = fit_a_tv.c1/2;
        subject_data.a_tv.mu    = fit_a_tv.b1;
        subject_data.a_tv.y     = fit_a_tv.a1*exp(-((xdata-fit_a_tv.b1).^2)/fit_a_tv.c1^2);
    catch err
        subject_data.a_tv.sigma = NaN;
        subject_data.a_tv.mu    = NaN;
        subject_data.a_tv.y     = zeros(13,1);        
    end

    try
        fit_t_av = fit(xdata.',subject_data.t_av.data.','gauss1');
        subject_data.t_av.sigma = fit_t_av.c1/2;
        subject_data.t_av.mu    = fit_t_av.b1;
        subject_data.t_av.y     = fit_t_av.a1*exp(-((xdata-fit_t_av.b1).^2)/fit_t_av.c1^2);
    catch err
        subject_data.t_av.sigma = NaN;
        subject_data.t_av.mu    = NaN;
        subject_data.t_av.y     = zeros(13,1);        
    end
    
    try
        fit_v_at = fit(xdata.',subject_data.v_at.data.','gauss1');
        subject_data.v_at.sigma = fit_v_at.c1/2;
        subject_data.v_at.mu    = fit_v_at.b1;
        subject_data.v_at.y     = fit_v_at.a1*exp(-((xdata-fit_v_at.b1).^2)/fit_v_at.c1^2);
    catch err
        subject_data.v_at.sigma = NaN;
        subject_data.v_at.mu    = NaN;
        subject_data.v_at.y     = zeros(13,1);        
    end    
end

function subject_data = setFit2(subject_data, xdata)
    [subject_data.a_tv.sigma, subject_data.a_tv.mu, subject_data.a_tv.y] = gFit(subject_data.a_tv.data, xdata);
    [subject_data.t_av.sigma, subject_data.t_av.mu, subject_data.t_av.y] = gFit(subject_data.t_av.data, xdata);
    [subject_data.v_at.sigma, subject_data.v_at.mu, subject_data.v_at.y] = gFit(subject_data.v_at.data, xdata);
end
