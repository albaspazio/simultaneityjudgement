addpath('bootstrap_fit'); 
addpath('classes'); 
addpath('utility'); 

%% edit this
root_dir            = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/';
experiment_folder   = 'blind_adult_iit_1'; ...'lv_children_chiossone'; ...'sighted_children_bolzaneto';
task_folder         = 'ATBSSU';
results_folder      = 'results';
result_postfix      = '';

%% init variables
xlabels         = {'-1200', '-800', '-400', '-300', '-200','-100','-50','0','50','100','200','300','400','800','1200'};
titleLabels     = {'A vs T'};
ylimits         = [0, 100];
xdata           = [-1200, -800, -400, -300, -200, -100, -50, 0, 50, 100, 200, 300, 400, 800, 1200];


%% start processing
data_dir        = fullfile(root_dir, experiment_folder, task_folder);

if ~exist(data_dir, 'dir')
    disp("ERROR.....DATA FOLDER DOES NOT EXIST!") 
    return;
end

result_file     = fullfile(pwd, results_folder, strcat(experiment_folder, "_", [task_folder result_postfix], ".dat"));

files           = dir(strcat(data_dir,'/*.txt'));
nsubj           = length(files);
subjects        = GroupATBSS(nsubj, xlabels, xdata, titleLabels, ylimits);

%% calculate the three timeseries (A_TV, T_AV, V_AT) of mean-yes-responses for each of the 13 delays.
... e.g. [0, 25, 25, 75, 100, 87.5,  91.6,  87.5, 87.5, 87.5, 62.5, 50,25]
for f=1:nsubj
    file_name       = files(f).name; ...'agga_8_1_ATVBSSU_TD_18122020_112553.txt';
    data            = tdfread(fullfile(data_dir, file_name));
    filename_parts  = split(file_name, '_');

    % create a Subject instance, parse and fit data
    subj            = SubjectATBSS( filename_parts{1}, ...              label
                                    str2double(filename_parts{2}), ...  age
                                    filename_parts{3}, ...              gender
                                    filename_parts{5}, ...              population
                                    filename_parts{4}, ...              experiment
                                    xdata, data);   ...                 data, xdata
    
    subjects        = subjects.add(subj);   % add subects to group
end

clear filename_parts file f files nsubj label

...subjects.plotSubject("ws", xdata, titleLabels)

...arni = subjects.getSubjectByLabel('ws');


...subjects.plotSubjectsSJ2("age", [8, 9]);

...subjects.plotSubjectsSJ2('label', 'ws');
subjects.plotSubjectsSJ2();
subjects.plotSubjectsGFit();
...subjects.plotSubjectsGFit("age", [8, 9]);
...subjects.plotSubjectsGFit("age", [8, 9], "gender", "m");

subjects.create_tabbed_data(result_file);
