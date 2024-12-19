function perceptualScore = PESDM_cal(input_hrir,ref_hrir,azimuth,mode)
% This function provides the calculation process for the Perceptually Enhanced Spectral Distance Metric for 
%        Head-Related Transfer Function Quality Prediction.
%
% SIMPLE USAGE EXAMPLES:
%   perceptualScore = PESDM_cal(input_hrir,ref_hrir,0,'L');
%
% INPUT:
%   input_hrir  =  HRIR to be evaluated, such as the HRIR obtained through modeling or prediction.
%   ref_hrir    =  reference HRIR.
%   azimuth     =  azimuth of the HRIR ranges from 0 to 360 degrees in a right-handed coordinate system, 
%                  where an azimuth of 90 degrees represents the subject's left side, and 270 degrees represents the right side.
%   mode        =  'L'/'R', the HRIR is for either the left ear ('L') or the right ear ('R').
%   NOTE: The sampling rate of the HRIRs needs to be 48,000 Hz.
% 
% OUTPUTS:
%   perceptualScore = perceptual difference between the input HRIR and the reference HRIR, 
%                     where a score of 100 indicates no difference, and a lower score indicates a greater difference.
%
% EXTERNAL DEPENDENCIES:
%   THE AUDITORY MODELING TOOLBOX 1.5.0 (https://amtoolbox.org/amt-1.5.0/doc/base/index.php)
%   Since the calculation of inner ear excitation patterns involves certain functions in AMT, it is necessary to load AMT.
% 
% REFERENCES:
%   Dingding Yao, Jiale Zhao, et al. "Perceptually Enhanced Spectral Distance Metric for Head-Related Transfer Function 
%        Quality Prediction." The Journal of the Acoustical Society of America.
%        
% AUTHOR: Dingding Yao - yaodingding(at)hccl.ioa.ac.cn
% December 2024

%% Generic parameters setting
% 0. Discrimination threshold: T 
JND_threshold=1.4; % The unit is dB (decibels).
% 1. Feature combination coefficient: alpha
alpha = 1.7;
% 2. Binaural weight parameter: Phi 
Phi   = 128;       % The unit is degrees.
% 3. Sensitivity coefficient: Gamma
Gamma = 0.0993;

%% HRIR information
Fs       = 48000; % Sampling rate (in Hz)
HRIR_len = 512;   % impulse response length

%% Spectral analysis
% MidEar_ir: The middle ear transfer function, derived from Lopez-Poveda, E. A., and Meddis, R. (2001).
MidEar_ir = middleearfilter(48000); % Url: http://amtoolbox.org/amt-1.5.0/doc/common/middleearfilter.php
ERB_len   = 1;      % The bandwidth of spectral analysis.
fre_flow  = 200;    % The lower frequency boundary of spectral analysis (in Hz).
fre_high  = 18000;  % The upper frequency boundary of spectral analysis (in Hz).
% Middle ear processing.
input = conv(input_hrir(getBeginPoint(input_hrir):end),MidEar_ir);
ref   = conv(ref_hrir(getBeginPoint(ref_hrir):end),MidEar_ir);        
input = input(:); ref = ref(:);
% Calculation of inner ear excitation patterns (EP).
[~, EP_input] = gammatone_fromAMT(input(1:HRIR_len),Fs,fre_flow,fre_high,ERB_len);
[~, EP_ref]   = gammatone_fromAMT(ref(1:HRIR_len),Fs,fre_flow,fre_high,ERB_len);
EP_input_db   = 20*log10(EP_input); % Excitation patterns of the input HRIR
EP_ref_db     = 20*log10(EP_ref); % Excitation patterns of the reference HRIR

%% Threshold discrimination
JND_index = find(abs(EP_input_db-EP_ref_db)<=JND_threshold); % Frequency band indices below the threshold.
EP_input_db(JND_index) = EP_ref_db(JND_index);

%% Feature combination
% d_MSE
d_MSE   = d_mse(EP_input_db,EP_ref_db);
% d_Std
d_Std   = d_std(EP_input_db,EP_ref_db);
% two distance metrics are combined to form the mixed distance
D_mixed = d_MSE+alpha*d_Std;

%% Binaural weighting
% azimuth needs to be converted to the range of -90 to 90 degrees.
% -90 degrees corresponds to the direction straight to the right, and 90 degrees corresponds to the direction straight to the left.
if azimuth>=0 && azimuth<=90
    azimuth_converted = azimuth;
end
% front-back mirroring
if azimuth>90 && azimuth<270
    azimuth_converted = 180-azimuth;
end
% azimuths from 270 degrees to 360 degrees are mapped to -90 degrees to 0 degrees.
if azimuth>=270 && azimuth<=360
    azimuth_converted = azimuth-360;
end
% determination of the left ear and right ear.
if strcmp(mode,'L')
    w_LR= (1+exp(-1*azimuth_converted/Phi)).^-1;
end
if strcmp(mode,'R')
    w_LR= 1-(1+exp(-1*azimuth_converted/Phi)).^-1;
end
% binaural weighting mixed distance
m_LR = w_LR.*(D_mixed);

%% Perceptual outcome estimation
perceptualScore = 100*exp(-1*Gamma*m_LR);

end

%% Time of arrival (TOA) detection.
function beginPoint = getBeginPoint(HRIR)
% Note: The threshold for leading-edge detection can range from 5% to 15%; this study used 10%. 
%       If the modeled or processed HRIR has poor stability, a threshold of 5% can be adopted to avoid energy loss.
threshold = 0.1*max(abs(HRIR)); % identifying the leading edge of the HRIR when it first reached 10% of its maximum peak amplitude
for k = 1 : length(HRIR)
    if abs(HRIR(k)) > threshold
        beginPoint = k;
        break;
    end
end
end

%% Inner ear excitation patterns calculation. Derived by rewriting BAUMGARTNER2014_SPECTRALANALYSIS
function [bands_fc, ep_out] = gammatone_fromAMT(sig,fs,flow,fhigh,space)
insig     = sig(:);
fc        = audspacebw(flow,fhigh,space,'erb');  % Url: http://ltfat.github.io/doc/auditory/audspacebw.html
[bgt,agt] = gammatone(fc,fs,'complex');          % Url: http://amtoolbox.org/amt-1.5.0/doc/common/gammatone.php
mp        = 2*real(ufilterbankz(bgt,agt,insig)); % Url: http://amtoolbox.org/amt-1.5.0/doc/common/ufilterbankz.php
% Averaging over time (RMS)
ep_out    = ((rms(mp)));  % rms out
bands_fc  = fc;
end

%% Calculation of the Mean Squared Error (MSE) between two input signals.
function out = d_mse(hrtf_1,hrtf_2)
% The units for hrtf_1 and hrtf_2 are dB.
hrtf_1  = hrtf_1(:);
hrtf_2  = hrtf_2(:);
freqlen = length(hrtf_1);
out     = sqrt((1/freqlen)*sum(power(hrtf_1-hrtf_2,2)));
end

%% Calculation of the standard deviation (Std) between two input signals.
function out = d_std(hrtf_1,hrtf_2)
% The units for hrtf_1 and hrtf_2 are dB.
hrtf_1 = hrtf_1(:);
hrtf_2 = hrtf_2(:);
out    = std(hrtf_1-hrtf_2);
end