clear all;
close all;
clc;
addpath(genpath('../../src'));


%-- Parameters
acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
phantom_type = 2;           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality
data_type = 1;              %-- 1 = IQ || 2 = RF
flag_display = 0;

%-- Parsing parameter choices
switch acquisition_type    
    case 1
        acquisition = 'simulation';
        acqui = 'simu';
        flag_simu = 1;
    case 2
        acquisition = 'experiments';
        acqui = 'expe';
        flag_simu = 0;
    otherwise       
        acquisition = 'simulation';
        acqui = 'simu';
        flag_simu = 1;
end
switch phantom_type    
    case 1
        phantom = 'resolution_distorsion';
    case 2
        phantom = 'contrast_speckle';
    otherwise       %-- Do deal with bad values
        phantom = 'resolution';
end
switch data_type    
    case 1
        data = 'iq';
    case 2
        data = 'rf';
    otherwise       %-- Do deal with bad values
        data = 'iq';        
end

%-- Create path to load corresponding files
path_dataset = ['../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_pht = ['../../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];

dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);
pht = us_phantom();
pht.read_file(path_pht);

flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);

x = {'Window-11 Uniform', 'Window-12 Uniform', 'Window-13 Uniform', 'Window-14 Uniform', ...
    'Window-15 Uniform', 'Window-16 Uniform', 'Window-17 Uniform', 'Window-18 Uniform', ...
    'Window-21 Uniform', 'Window-22 Uniform', 'Window-31 Uniform', 'Window-32 Uniform', ...
    'Window-33 Uniform', 'Window-34 Uniform', 'Window-41 Uniform', 'Window-42 Uniform'};



for k=x
    s1 = ['Results weighting strategy 1/simulation/NO-MEAN-FILTER/',k{1}, '/TwoMinus/windowed_15.hdf5'];
    s1_txt = ['Results weighting strategy 1/simulation/NO-MEAN-FILTER/',k{1}, '/TwoMinus/windowed_15.txt'];

    s2 = ['Results weighting strategy 2/simulation/NO-MEAN-FILTER/',k{1},'/TwoMinus/windowed_15.hdf5'];
    s2_txt = ['Results weighting strategy 2/simulation/NO-MEAN-FILTER/',k{1}, '/TwoMinus/windowed_15.txt'];
    
    image1 = us_image();
    image1.read_file(s1);
    con1 = tools.contrast_score(scan, pht, image1, 9);
    
    tools.exec_evaluation_contrast_speckle(path_scan,path_pht, ...
        s1,flag_simu,flag_display,s1_txt);

    tools.exec_evaluation_contrast_speckle(path_scan,path_pht, ...
        s2,flag_simu,flag_display,s2_txt);
    
    

    image2 = us_image();
    image2.read_file(s2);
    con2 = tools.contrast_score(scan, pht, image2, 9);

    disp(mean(con1));
    disp(mean(con2));
end
