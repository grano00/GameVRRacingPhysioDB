function [fiducial_points, BL]=embdac(ECGLEAD,SAMPLING_RATE,MICROV_FOR_ADCUNIT)
% fiducial_points=embdac(ECGLEAD,SAMPLING_RATE,MICROV_FOR_ADCUNIT)
% Uses a slightly (very slightly) modified versions of OSEA
% http://www.eplimited.com
% to detect classical ECG fiducial points.
%
% ECGLEAD            : a single ECG lead (int16)
% SAMPLING_RATE      : the sampling rate in Hz
% MICROV_FOR_ADCUNIT : the number of microvolt which correspond to a single 
%                      ADC unit (floating point number)
%
% Each line in fiducial_points refers to a single beat
% fiducial_points(:,1) : beat fiducial point
% fiducial_points(:,2) : beat label
% fiducial_points(:,3) : qrs onset
% fiducial_points(:,4) : qrs offset
% fiducial_points(:,5) : beat onset
% fiducial_points(:,6) : beat onset
%
% beat label is one of the followings:
% (int)'N' : normal beat
% (int)'V' : ventricular ectopic beat
% (int)'A' : artifact

if( exist('embdac_core','file') ~= 3 ),
    if( (exist('embdac_core.c','file')~=2) || ...
            (exist('qrsdet.h','file')~=2) || ...
            (exist('qrsdet.c','file')~=2) || ... % era qrsdet2.c in precedenza ...
            (exist('qrsfilt.c','file')~=2) || ...
            (exist('bdac.c','file')~=2) || ...
            (exist('classify.c','file')~=2) || ...
            (exist('rythmchk.h','file')~=2) || ...
            (exist('analbeat.c','file')~=2) || ...
            (exist('match.c','file')~=2) || ...
            (exist('postclas.c','file')~=2) || ...
            (exist('bdac.h','file')~=2) || ...
            (exist('ecgcodes.h','file')~=2) || ...
            (exist('match.h','file')~=2) || ...
            (exist('rythmchk.h','file')~=2) || ...
            (exist('analbeat.h','file')~=2) || ...
            (exist('postclas.h','file')~=2) ),
        error('At least a file is missing.');
    else
        % mex embdac_core.c qrsdet2.c qrsfilt.c bdac.c classify.c rythmchk.c noisechk.c analbeat.c match.c postclas.c
        mex embdac_core.c qrsdet.c qrsfilt.c bdac.c classify.c rythmchk.c noisechk.c analbeat.c match.c postclas.c
    end,
end,

% The mex file does exist.

% OSEA expects the signal to be int16 and have microv_adcunit=5
% When it is 10uV, is ecgSample=raw_ecg*2 correct? YES
% When it is 1uV, is ecgSample=raw_ecg/5 correct? YES
mlead=double(ECGLEAD)*MICROV_FOR_ADCUNIT/5;

% OSEA expects the signal to be sampled at 200 Hz
[n_rat,d_rat]=rat(200/SAMPLING_RATE);
mlead=resample(mlead,n_rat,d_rat);
mlead=int16(round(mlead));

%% Call OSEA
BL=embdac_core(mlead,5,200);

%% The output is remapped to the original sampling rate
% the random rounding make the error zero on average, but it is an
% unneccessary computational burden...
% fiducial_points=round((BL-1)/n_rat*d_rat+1+1e-5*(rand(size(BL))-0.5));
fiducial_points=round((BL-1)/n_rat*d_rat+1);
fiducial_points(fiducial_points(:)<1)=nan;
fiducial_points(fiducial_points(:)>length(ECGLEAD))=nan;
fiducial_points(:,2)=BL(:,2);

return

