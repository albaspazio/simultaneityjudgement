%% edit this

root_dir = '/data/Dropbox/BIODOCS/projects/Apps/PsySuite/DATA/';
experiment_folder   = 'bambini_bolzaneto';   ... 'bambini_bolzaneto' 'adulti_iit_1'
task_folder = 'ATVBSSB';  ... 'ATBSSU' 'ATVBDSU' 'ATVBSSU'


new_header          = 'id\ttype\tdummy\tdelay\tanswer\tsuccess\telapsed';

%% start processing
data_dir    = fullfile(root_dir, experiment_folder, task_folder);
cd(data_dir)
files = dir('*.txt');
nsubj = length(files);

for f=1:nsubj

    file            = files(f).name;
        
    fid = fopen(file,'r');        % Open File to read
    i = 1;
    tline = 's';
    A = {[]};
    while ischar(tline)
        tline = fgetl(fid);
        
        if i>1
            A{i}=tline;
        end
        i = i+1;
    end
    
    fclose(fid);
    fid2=fopen(file,'w');            % Open file to write
    
    firstline = sprintf(new_header);
    fprintf(fid2,['%s',char([13,10])],firstline);
    for i=2:length(A)-1
        fprintf(fid2,['%s',char([13,10])],A{i});
    end
    
    fclose(fid2);
end
