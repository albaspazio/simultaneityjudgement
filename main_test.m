addpath('bootstrap_fit'); 
addpath('classes'); 
addpath('utility'); 

%% edit this
root_dir            = fullfile(pwd, 'input_data');
experiment_folder   = '';           ... possible subfolder of 'input_data'
task_folder         = 'ATVBSSU';    ... possible subfolder of experiment_folder
result_postfix      = '';           ... string to append to result file

%% group variables
xlabels         = {'-1200', '-800', '-400', '-300', '-200','-100','-50','0','50','100','200','300','400','800','1200'};
xdata           = [-1200, -800, -400, -300, -200, -100, -50, 0, 50, 100, 200, 300, 400, 800, 1200];
titleLabels     = {'A vs TV', 'T vs AV', 'V vs AT'};
ylimits         = [0, 100];

%% start processing

data_dir        = fullfile(root_dir, experiment_folder, task_folder);

if ~exist(data_dir, 'dir')
    disp("ERROR.....DATA FOLDER DOES NOT EXIST!") 
    return;
end

result_file     = fullfile(pwd, results_folder, strcat(experiment_folder, "_", [task_folder result_postfix], ".dat"));
errors_file     = fullfile(pwd, results_folder, strcat(experiment_folder, "_", [task_folder '_errors' result_postfix], ".dat"));

files           = dir(strcat(data_dir,'/*.txt'));   % list files contained within given folder
nsubj           = length(files);

subjects        = GroupATVB(nsubj, xlabels, xdata, titleLabels, ylimits); % class managing the subjects list

%% calculate the 1-3 timeseries of yes-responses for each of the 13/15 delays.
... e.g. [0, 25, 25, 75, 100, 87.5,  91.6,  87.5, 87.5, 87.5, 62.5, 50,25]
... input file names must be formatted as:  subjlabel_age_gender_experiment_population_date_time.txt
... e.g: 'agga_8_m_ATVBSSU_TD_18122020_112553.txt';

for f=1:nsubj
    file_name       = files(f).name; ...
    data            = tdfread(fullfile(data_dir, file_name));
    filename_parts  = split(file_name, '_');

    % create a Subject instance, parse and fit data
    subj            = SubjectATVB(  filename_parts{1}, ...              subjlabel
                                    str2double(filename_parts{2}), ...  age
                                    filename_parts{3}, ...              gender
                                    filename_parts{5}, ...              population
                                    filename_parts{4}, ...              experiment
                                    xdata, data);   ...                 xdata, data
    
    subjects        = subjects.add(subj);   % add subects to group

end

%% now data are loaded in subjects, an instance of GroupATVB, you can call Group's methods

...subjects.plotSubject("arni", xdata, titleLabels)         ... plot gaussian fit of a single subect by its label
...subjects.plotSubjectsSJ2("age", [8, 9]);                 ... plot and SJ2 (stanley) metrics of all subject aged 8 & 9
...subjects.plotSubjectsSJ2('label', 'arni');               ... plot and SJ2 (stanley) metrics of one subject
...subjects.plotSubjectsSJ2();                              ... plot and SJ2 (stanley) metrics of all subjects

...subjects.plotSubjectsGFit();                             ... plot gaussian fit of all subjects
...subjects.plotSubjectsGFit("age", [8, 9]);                ... plot gaussian fit of all subject aged 8 & 9
...subjects.plotSubjectsGFit("age", [8, 9], "gender", "m"); ... plot gaussian fit of all male subject aged 8 & 9

%% or you can get a specific subject instance and work on it
...arni = subjects.getSubjectByLabel('arni');           ... how to get a specific subject instance


%% create results file
subjects.create_tabbed_data(result_file, errors_file);               ... create results file

%% clear unnecessary variables
clear f file_name data filename_parts subj files
clear xlabels xdata ylimits titleLabels taskfolder root_dir data_dir results_folder result_file experiment_folder task_folder nsubj
