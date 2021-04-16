
%% SETUP THE MATLAB PATHS AND FILE NAMES
% make sure that fieldtrip and spm are not in your matlab path

mydir = '/Users/mortenk/Documents/MATLAB/'; % adapt to yours
workingdir = [mydir 'preprocess_HCP/'];
codedir = [mydir 'osl/'];
osldir = [codedir 'osl-core/'];
hmmdir = [codedir 'HMM-MAR/']; % latest version on git
datadrive = [mydir 'HCP/'];
%datadirRAW = [datadrive '/datasets/HCP_RAW/HCP_DATA/'];
datadirRECONS = [datadrive];
datadirOUT = [datadrive 'ts/HCP_MEG/fMRI_recons_perband/'];



%session_names = {'Restin','Motort','StoryM','Wrkmem'};
session_dirs = {'REST/'};   
%session_dirs = {'REST/','TASKS/LANGUAGE/',...
%    'TASKS/MOTOR/','TASKS/WORKMEM/'};   
%band_names = {'Slow','Alpha','Beta','Gamma'};
band_names = {'Broad'};

%orthogonalisation = 'innovations_mar';
%innovations_mar_order = 14;

addpath(osldir)
osl_startup

addpath([codedir '/ohba-external/ohba_utils'])
addpath([codedir '/ohba-external/nifti_tools'])
addpath(genpath([codedir '/spm12']))
addpath(genpath(hmmdir))

parcdir = [codedir 'parcellations/'];
%parcfile = [parcdir '/Parcels_Combo_8mm_OE.nii.gz'];
%parcfile = [parcdir 'dk_cortical.nii.gz'];
parcfile = [parcdir 'dkt62_cortical_4D_8mm.nii.gz'];
%parcfile = [parcdir 'dbs80symm_8mm_4D.nii.gz'];
%parcfile = [parcdir 'fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz'];

maskfname = [codedir 'std_masks/MNI152_T1_8mm_brain_mask.nii.gz'];
%p_gordon = parcellation(parcfile);

parcprefix = 'dbs80_';
% PCA_weights = dlmread([workingdir 'PCA_comp.txt']);

% p_std = parcellation(fullfile(codedir,'parcellations',...
%     'fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'));


%LANGDIR = '/data/analysis_data/datasets/HCP/TASKS/LANGUAGE/';
%dat = spm_eeg_load([LANGDIR '158136_MEG_9-StoryM_tmegpreproc_BSENT']);

% %% Preproc
% 
% Fs = 200;
% Bands = [0 7; 7.5 13; 14 30; 31 Inf];
% 
% Done = cell(4,1);
% 
% % Cycle through fif files
% for j = 2
%     switch j
%         case 1, datadir = [datadirRECONS session_dirs{j} '/*dat' ];
%         case 2, datadir = [datadirRECONS session_dirs{j} '/*BU.dat' ];
%         case 3, datadir = [datadirRECONS session_dirs{j} '/*.dat' ];
%         case 4, datadir = [datadirRECONS session_dirs{j} '/*TIM.dat' ];
%     end
%     files = dir(datadir);
%     datadirOUT_j = [datadirOUT session_dirs{j}];
%     Done{j} = false(length(files),1);
%     for s = 1:length(files)
%         file_in = [datadirRECONS session_dirs{j} files(s).name];
%         EXISTS = true;
%         for i = 1:size(Bands,1)
%             file_out = [datadirOUT_j band_names{i} '/' files(s).name];
%             EXISTS = EXISTS & exist(file_out);
%         end
%         if EXISTS, Done{j}(s) = true; continue; end
%         try
%             % Read data
%             D = spm_eeg_load(file_in);
%             D = D.montage('switch',2);
%             D = osl_resample_meeg(D,'ft_8mm_brain_mask.nii.gz','MNI152_T1_8mm_brain_mask');
%             D = D.montage('remove',2);
% %             for j1=1:size(D,3)
% %                 for j2=1:size(D,2)
% %                     which_nan = isnan(D(:,j2,j1));
% %                     D(which_nan,j2,j1) = 0; 
% %                 end
% %             end
%             % parcellate
%             S                   = [];
%             S.D                 = D;
%             S.parcellation      = parcfile;
%             S.orthogonalisation = 'none';
%             %S.innovations_mar_order = innovations_mar_order;
%             S.whichmethod            = 'spatialBasis';
%             S.normalise_voxeldata = 0;
%             S.prefix = parcprefix;
%             S.maskfname = maskfname;
%             [p_D,parcelWeights,parcelAssignments] = osl_apply_parcellation(S);
%             p_D.parcellation.weights = parcelWeights;
%             p_D.parcellation.assignments = parcelAssignments;
%             %p_D.save;
%             Xraw = spm_eeg_load(p_D);
%             Xraw = Xraw(:,:,:);
%             Xraw = reshape(Xraw,[p_D.nchannels,p_D.nsamples*p_D.ntrials])';
%             T = p_D.nsamples * ones(p_D.ntrials,1);
%             %[Xraw,T] = read_spm_file(p_D);
%             % Filter into three bands
%             for i = 1:size(Xraw,2)
%                 are_nan = isnan(Xraw(:,i));
%                 Xraw(are_nan,i) = mean(Xraw(~are_nan,i));
%             end
%             non_zero = var(Xraw)>0;
%             Xraw = Xraw(:,non_zero);
%             for i = 1:size(Bands,1)
%                 file_out = [datadirOUT_j band_names{i} '/' files(s).name];
%                 X = filterdata(Xraw,T,Fs,Bands(i,:));
%                 X = single(rawsignal2power(X,T) * PCA_weights(non_zero,4:8));
%                 save(file_out,'X','T')
%             end
%             Done{j}(s) = true;
%         catch
%             warning(['Something went wrong with task ' ...
%                 num2str(j) ' session ', num2str(s)])
%         end
%         disp(['Done task ' num2str(j), ' subject ' num2str(s)])
%     end
%     system(['rm ' datadirRECONS session_dirs{j} 'Gordon_' files(s).name '*']) 
% end
% 
% save([workingdir '/Done.mat'],'Done')
% 

%% Do the same but without the PCA reduction

Fs = 200;
%Bands = [0 7; 7.5 13; 14 30; 31 Inf];
Bands = [1 98];

Done_2 = cell(4,1);

% Cycle through fif files
for j = [1 4]
    switch j
        case 1, datadir = [datadirRECONS session_dirs{j} '/*mat' ];
        case 2, datadir = [datadirRECONS session_dirs{j} '/*BU.dat' ];
        case 3, datadir = [datadirRECONS session_dirs{j} '/*.dat' ];
        case 4, datadir = [datadirRECONS session_dirs{j} '/*TIM.dat' ];
    end
    files = dir(datadir);
    datadirOUT_j = [datadirOUT session_dirs{j}];
    Done_2{j} = false(length(files),1);
    for s = 1:length(files)
        file_in = [datadirRECONS session_dirs{j} files(s).name];
        EXISTS = true;
        for i = 1:size(Bands,1)
            file_out = [datadirOUT_j band_names{i} '/' files(s).name];
            EXISTS = EXISTS & exist(file_out);
        end
        if EXISTS, Done_2{j}(s) = true; continue; end
        try
            % Read data
            if strcmp(files(s).name(1:6),'dbs80'), continue; end % residual from previous run
            disp(['Task ' num2str(j), ' subject ' num2str(s)])
            D = spm_eeg_load(file_in);
            D = D.montage('switch',2);
            D = osl_resample_meeg(D,'ft_8mm_brain_mask.nii.gz','MNI152_T1_8mm_brain_mask.nii.gz');
            D = D.montage('remove',2);
            %             for j1=1:size(D,3)
            %                 for j2=1:size(D,2)
            %                     which_nan = isnan(D(:,j2,j1));
            %                     D(which_nan,j2,j1) = 0;
            %                 end
            %             end
            % parcellate
            S                   = [];
            S.D                 = D;
            S.parcellation      = parcfile;
            S.orthogonalisation = 'none';
            %S.innovations_mar_order = innovations_mar_order;
            S.whichmethod            = 'spatialBasis';
            S.normalise_voxeldata = 0;
            S.prefix = parcprefix;
            S.maskfname = maskfname;
            [p_D,parcelWeights,parcelAssignments] = osl_apply_parcellation(S);
            p_D.parcellation.weights = parcelWeights;
            p_D.parcellation.assignments = parcelAssignments;
            %p_D.save;
            Xraw = spm_eeg_load(p_D);
            Xraw = Xraw(:,:,:);
            Xraw = reshape(Xraw,[p_D.nchannels,p_D.nsamples*p_D.ntrials])';
            T = p_D.nsamples * ones(p_D.ntrials,1);
            %[Xraw,T] = read_spm_file(p_D);
            % Filter into three bands
            for i = 1:size(Xraw,2)
                are_nan = isnan(Xraw(:,i));
                Xraw(are_nan,i) = mean(Xraw(~are_nan,i));
            end
            non_zero = var(Xraw)>0;
            Xraw = Xraw(:,non_zero);
            for i = 1:size(Bands,1)
                file_out = [datadirOUT_j band_names{i} '/' files(s).name];
                X = filterdata(Xraw,T,Fs,Bands(i,:));
                X = single(rawsignal2power(X,T));
                %X = single(rawsignal2power(X,T) * PCA_weights(non_zero,4:8));
                save(file_out,'X','T')
            end
            Done_2{j}(s) = true;
        catch
            warning(['Something went wrong with task ' ...
                num2str(j) ' session ', num2str(s)])
        end
        disp(['Done task ' num2str(j), ' subject ' num2str(s)])
    end
    %system(['rm ' datadirRECONS session_dirs{j} 'Gordon_' files(s).name '*'])
end

save([workingdir '/Done.mat'],'Done_2','--append')
    
    
    
    

