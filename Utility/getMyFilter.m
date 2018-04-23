function [ f ] = getMyFilter( type,  passband, stopband, sr, attenuation, method)
%GETFILTER this function return the differents types of FIR filter.
%   Input argument are: 
%       - type: 'lp', 'hp', 'bp', 'bs'
%       - passband and stopband could be scalar values for lp and hp, or vectors of two
%       elements for bp and bs. They define, respectively, the starting
%       value and the stopping value of the filter cutoff.
%       - sr is the sample rate of the data
%       - attenuation: scalar in db that express the stopband attuation
%       - method: the filtering method. The standard method is equiripple.
%       The other methods are:
%           - 'cls' designs an FIR filter using constrained least squares. The method minimizes the discrepancy between a specified arbitrary piecewise-linear function and the filter’s magnitude response. At the same time, it lets you set constraints on the passband ripple and stopband attenuation.
%           - 'equiripple' designs an equiripple FIR filter using the Parks-McClellan algorithm. Equiripple filters have a frequency response that minimizes the maximum ripple magnitude over all bands.
%           - 'freqsamp' designs an FIR filter of arbitrary magnitude response by sampling the frequency response uniformly and taking the inverse Fourier transform.
%           - 'kaiserwin' designs an FIR filter using the Kaiser window method. The method truncates the impulse response of an ideal filter and uses a Kaiser window to attenuate the resulting truncation oscillations.
%           - 'ls' designs an FIR filter using least squares. The method minimizes the discrepancy between a specified arbitrary piecewise-linear function and the filter’s magnitude response.
%           - 'maxflat' designs a maximally flat FIR filter. These filters have a smooth monotonic frequency response that is maximally flat in the passband.
%           - 'window' uses a least-squares approximation to compute the filter coefficients and then smooths the impulse response with 'Window'.

f = -1;

%define the correct names for the type of filter
switch type
    case 'lp'
        type = 'lowpassfir';
    case 'hp'
        type = 'highpassfir';
    case 'bp'
        type = 'bandpassfir';
    case 'bs'
        type = 'bandstopfir';
    otherwise
        error('Insert a valid type of filter');
end

%check if passband and stopband are correctly setted
if((strcmp(type,'lowpassfir') || strcmp(type,'highpassfir')) && (~isscalar(passband) || ~isscalar(stopband)))
    error('insert a scalar value for lp and hp in passband and stopband');
end

if((strcmp(type,'bandpassfir') || strcmp(type,'bandstopfir')) && (length(passband) ~= 2 || length(stopband) ~= 2))
    error('insert a 2d vector the value for bp and bs in passband and stopband');
end


%
if (~exist('method','var'))
    method = 'equiripple';
end
%TOCOMPLETE THE OTHERS CHECKS

if((strcmp(type,'lowpassfir') || strcmp(type,'highpassfir')))
   f = designfilt(type, 'PassbandFrequency', passband, 'StopbandFrequency', stopband, ...
        'PassbandRipple', 1, 'StopbandAttenuation', attenuation, ...
        'SampleRate', sr, 'DesignMethod', method);
end

if((strcmp(type,'bandpassfir') || strcmp(type,'bandstopfir')))
    f = designfilt(type, 'PassbandFrequency1', passband(1), 'PassbandFrequency2', passband(2), ...
        'StopbandFrequency1', stopband(1), 'StopbandFrequency2', stopband(2), ...
        'StopbandAttenuation1', attenuation, 'PassbandRipple', 1, ...
        'StopbandAttenuation2', attenuation, 'SampleRate', sr, 'DesignMethod', method);
end
end

