%% read the json associated to each result file, extract age and gender and append to results' filename
% from  : hega_ATVBSSU_TD_15122020_084435.txt
% to    : hega_age_gender_ATVBSSU_TD_15122020_084435.txt

%% edit this
% root_dir    = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/bambini_bolzaneto/';
% root_dir    = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/nonvedenti_ragazzi_chiossone/';
root_dir    = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/blind_adult_iit_1/';
...data_dir        = 'E:\\Dropbox\\BIODOCS\\projects\\Apps\\PsySuite\\DATA\\bambini_bolzaneto\\';
    

task_folder = 'ATBSSU'; ...'ATVBDSU';  ... 'ATBSSU' 'ATVBDSU' 'ATVBSSU'


data_dir    = fullfile(root_dir, task_folder);
json_dir    = fullfile(data_dir, 'json');
%% start processing

curr_dir = pwd;
% read results
cd(data_dir);
res_files = dir('*.txt');
nsubj = length(res_files);

% read json
cd(json_dir);
json_files = dir('*.json');
njsons = length(json_files);

cd(curr_dir);
try
    for f=1:nsubj

        res_file    = res_files(f).name;
        fid         = fopen(fullfile(data_dir, res_file),'r');        % Open File to read

        arr_subjname    = split(res_file, "_");
        subjname        = arr_subjname{1};
        
        for s=1:njsons
            
            % get json file first part (subject label)
            arr_jsonname    = split(json_files(s).name, "_");
            jsonname        = arr_subjname{1};

            if strcmp(subjname, jsonname) == 1

                [age, gender] = get_age_gender(fullfile(json_dir, json_files(s).name), subjname);
                if strcmp(age, "") == 0

                    [~, f,ext] = fileparts(res_file);
                    new_name = strcat(subjname,'_', age, '_', gender, '_', strjoin(arr_subjname(2:end),'_')) ; 
                    movefile(fullfile(data_dir, res_file), fullfile(data_dir, new_name));                 
                    break;
                end
            end
        end

        fclose(fid);
    end
catch err
    err
end

function [age, gender] = get_age_gender(filename, subjname)
    age     = '';
    gender  = '';
    
    fid = fopen(filename); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);    
    if strcmp(val.label, subjname)
        age     = num2str(val.age);
        
        if val.gender == 0
            gender = 'm';
        else
            gender  = 'f';
        end
    end
end