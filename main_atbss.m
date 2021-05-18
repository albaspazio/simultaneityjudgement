addpath('bootstrap_fit'); 
addpath('classes'); 
addpath('utility'); 

%% edit this
root_dir            = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/';
experiment_folder   = 'bambini_bolzaneto';   ... 'bambini_bolzaneto' 'adulti_iit_1'
task_folder         = 'ATBSS';
results_folder      = 'results';


%% init variables
xlabels         = {'-1200', '-800', '-400', '-300', '-200','-100','-50','0','50','100','200','300','400','800','1200'};
titleLabels     = {'A vs T'};
ylimits         = [0, 100];
xdata           = [-1200, -800, -400, -300, -200, -100, -50, 0, 50, 100, 200, 300, 400, 800, 1200];


%% start processing
data_dir        = fullfile(root_dir, experiment_folder, task_folder);
result_file     = fullfile(pwd, results_folder, strcat(experiment_folder, "_", task_folder, ".dat"));

files           = dir(strcat(data_dir,'/*.txt'));
nsubj           = length(files);
subjects        = GroupATBSS(nsubj, xlabels, xdata, titleLabels, ylimits);

%% calculate the three timeseries (A_TV, T_AV, V_AT) of mean-yes-responses for each of the 13 delays.
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

...subjects.plotSubject("arni", xdata, titleLabels)
subjects.plotSubjects();

subjects.create_tabbed_data(result_file);
