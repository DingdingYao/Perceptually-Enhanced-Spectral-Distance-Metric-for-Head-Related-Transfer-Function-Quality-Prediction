% This script performs HRTF quality prediction example. The data used is from Dingding Yao, Jiale Zhao, et al.
%   "Perceptually Enhanced Spectral Distance Metric for Head-Related Transfer Function Quality Prediction,"
%   The Journal of the Acoustical Society of America.
%   It includes data from 16 subjects, with each dataset containing distorted HRIRs
%   across 32 spatial directions using 7 common HRTF processing methods, with 6 varying levels of LSD,
%   as well as the corresponding MUSHRA localization and timbre performance of these distorted HRIRs.
%
% EXTERNAL DEPENDENCIES:
%   THE AUDITORY MODELING TOOLBOX 1.5.0 (https://amtoolbox.org/amt-1.5.0/doc/base/index.php)
%
% AUTHOR: Dingding Yao - yaodingding(at)hccl.ioa.ac.cn
% December 2024

close all;clc;clear;
amt_start; % the calculation of inner ear excitation patterns involves certain functions in AMT

%% Subject ID
%   Subject ID description:
%   'L' indicates the subjects only distorted the HRTF of their left ear while keeping the HRTF of the right ear unchanged.
%   'R' indicates the subjects did the opposite.
%    8 subjects from 'L' group: 'L1', 'L2','L3','L4','L5','L6','L7','L8';
%    8 subjects from 'R' group: 'R1', 'R2','R3','R4','R5','R6','R7','R8'.
subject_IDs = {'L1'};% 'L1', 'L2','L3','L4','L5','L6','L7','L8','R1', 'R2','R3','R4','R5','R6','R7','R8'

%% HRTF quality prediction
for id = 1:length(subject_IDs)
    %% Loading the listening test data and results.
    %   The file x.mat contains the following data:
    %   1. HRIRs_tested         (32 dirs x 7 methods x [1 hid. ref. + 1 hid. anc. + 6 LSDs] x 2 ears x 512 irlen)
    %      -- Distorted HRIRs under corresponding experimental conditions.
    %   2. Results_localization (32 dirs x 7 methods x [1 hid. ref. + 1 hid. anc. + 6 LSDs])
    %      -- Similarity scores in perceived location between distorted HRIRs and reference HRIRs.
    %   3. Results_timbre       (32 dirs x 7 methods x [1 hid. ref. + 1 hid. anc. + 6 LSDs])
    %      -- Similarity scores in perceived coloration between distorted HRIRs and reference HRIRs.
    %   4. Dirs                 (32 dirs x [azimuth elevation])
    %   5. Methods_list         (7 HRTF processing methods for distortion)
    %      -- 'FSS', 'FSC', 'ARMA', 'PCA', 'SH', 'SELECTION', 'DL'.
    %   6. Fs: sampling rate of HRIRs.
    %   7. Mode_LR: 'L'/'R', the reference HRIR is for either the left ear ('L') or the right ear ('R').
    load(sprintf('Result_data\\%s.mat',subject_IDs{id}),'HRIRs_tested','Results_localization','Results_timbre','Dirs','Methods_list','Fs','Mode_LR');
    if strcmp(Mode_LR,'L')
        ear=1;
    end
    if strcmp(Mode_LR,'R')
        ear=2;
    end
    
    %% Extract valid data.
    Results_everyDir = [];
    for i = 1:size(Dirs) % 32 dirs
        Results_everyDir{i}.averaged_results = []; % Variable used to store the averaged perceptual results of localization and timbre.
        Results_everyDir{i}.distorted_hrirs  = []; % Variable used to store the corresponding distorted HRIRs.
        Results_everyDir{i}.target_hrir      = squeeze(HRIRs_tested(i,1,1,ear,:)); % Extract the target HRIR, i.e., the reference HRIR.
        counter = 1;
        for j = 1:length(Methods_list) % 7 HRTF processing methods
            ScoreThreshold = 90; % The trials with rating scores greater than 90 on the hidden reference were considered valid data.
            if (Results_localization(i,j,1))>ScoreThreshold && (Results_timbre(i,j,1))>ScoreThreshold
                for k = 3:8 % [k=3:8] corresponds to LSD levels [lsd=1:6]; k=1 corresponds to the hid. ref.; k=2 corresponds to the hid. anc.
                    Results_everyDir{i}.averaged_results(counter)  = (Results_localization(i,j,k)+Results_timbre(i,j,k))/2;
                    Results_everyDir{i}.distorted_hrirs(counter,:) = squeeze(HRIRs_tested(i,j,k,ear,:));
                    counter = counter + 1;
                end
            end
        end
    end
    
    %% Test PESDM
    Score_predicted = []; % Variable used to store the predicted perceptual results.
    Score_measured  = []; % Variable used to store the measured perceptual results.
    for i = 1:length(Results_everyDir)
        % load valid data
        target_hrir         = Results_everyDir{1, i}.target_hrir;
        distorted_hrirs     = Results_everyDir{1, i}.distorted_hrirs;
        Score_measured_temp = Results_everyDir{1, i}.averaged_results;
        for j = 1:length(Score_measured_temp)
            Score_predicted_temp(j) = PESDM_cal(distorted_hrirs(j,:),target_hrir,Dirs(i,1),Mode_LR);
        end
        % data combine
        Score_predicted = [Score_predicted; Score_predicted_temp'];
        Score_measured  = [Score_measured;  Score_measured_temp'];
        % free temporary data
        Score_predicted_temp = [];
        Score_measured_temp = [];
    end
    % Calculate the correlation coefficient.
    corr_proposed(id) = corr(Score_measured,Score_predicted);
    % Plot the relationships between the results predicted by PESDM and the measured results.
    figure;
    plot(Score_measured,100-Score_predicted,'.k','markersize', 14);
    lm = fitlm(Score_measured,100-Score_predicted);
    hold on;plot(Score_measured, predict(lm, Score_measured), 'r','Linewidth', 2);
    ylim([0 75]);
    xlim([25 100]);
    set(gca,'XTick',25:25:100);
    set(gca, 'YTick', 0:25:75);
    set(gca, 'YTickLabel',{'100','75','50','25'});
    ax = gca;
    set(ax, 'FontSize', 16);
    xlabel('Perceptual results','fontsize', 18);
    ylabel('Predicted Scores','fontsize', 18);
    title(sprintf('%s, |r| = %.2f','Proposed method',abs(corr_proposed(id))),'fontsize', 22);pbaspect([1.28 1 1])
    
end

