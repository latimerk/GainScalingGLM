function [cwtx] = wcoherence_getWaveletTransform(x,varargin)
%Wavelet coherence
% WCOH = WCOHERENCE(X,Y) returns the magnitude-squared wavelet coherence
% between the equal-length 1-D real-valued signals X and Y using the
% analytic Morlet wavelet. X and Y must have at least 4 samples. The
% wavelet coherence is computed over logarithmic scales using 12 voices per
% octave. The number of octaves is equal to floor(log2(numel(X)))-1.
%
% [WCOH,WCS] = WCOHERENCE(X,Y) returns the wavelet cross spectrum of X and
% Y in WCS.
%
% [WCOH,WCS,PERIOD] = WCOHERENCE(X,Y,Ts) uses the positive <a href="matlab:help duration">duration</a>, Ts,
% to compute the scale-to-period conversion, PERIOD. PERIOD is
% an array of durations with the same Format property as Ts.
%
% [WCOH,WCS,F] = WCOHERENCE(X,Y,Fs) uses the positive sampling frequency,
% Fs, in hertz to compute the scale-to-frequency conversion, F. If you
% output F without specifying a sampling frequency, WCOHERENCE uses
% normalized frequency in cycles/sample. The Nyquist frequency is 1/2.
% You cannot specify both a sampling frequency and a duration.
%
% [WCOH,WCS,F,COI] = WCOHERENCE(...) returns the cone of influence in
% cycles/sample for the wavelet coherence. If you specify a sampling
% frequency, Fs, in hertz, the cone of influence is returned in hertz.
%
% [WCOH,WCS,PERIOD,COI] = WCOHERENCE(...,Ts) returns the cone of influence
% in periods for the wavelet coherence. Ts is a positive <a href="matlab:help duration">duration</a>. COI is
% an array of durations with same Format property as Ts.
%
% [...] = WCOHERENCE(...,'VoicesPerOctave',NV) specifies the number of
% voices per octave to use in the wavelet coherence. NV is an integer in
% the range [10,32].
%
% [...] = WCOHERENCE(...,'NumScalesToSmooth', NS) specifies the number of
% scales to smooth as a positive integer less than one half the number of
% scales. If unspecified, NS defaults to the number of voices per octave. A
% moving average filter is used to smooth across scale.
%
% [...] = WCOHERENCE(...,'NumOctaves',NOCT) specifies the number of
% octaves to use in the wavelet coherence. NOCT is a positive integer
% between 1 and floor(log2(numel(X)))-1. If unspecified, NOCT defaults to
% floor(log2(numel(X)))-1.
%
% WCOHERENCE(...) with no output arguments plots the wavelet coherence in
% the current figure window along with the cone of influence. For areas
% where the coherence exceeds 0.5, arrows are also plotted to show the
% phase lag between X and Y. The phase is plotted as the lag between Y and
% X. The arrows are spaced in time and scale. 
%
% WCOHERENCE(...,'PhaseDisplayThreshold',PT) displays phase vectors for
% regions of coherence greater than or equal to PT. PT is a real-valued
% scalar between 0 and 1. This name-value pair is ignored if you call
% WCOHERENCE with output arguments.
%
%   % Example 1:
%   %   Plot the wavelet coherence for two signals. Both signals consist
%   %   of two sine waves (10 and 50 Hz) in white noise. The sine waves
%   %   have different time supports. The sampling interval frequency is
%   %   1000 Hz.
%   t = 0:0.001:2;
%   x = cos(2*pi*10*t).*(t>=0.5 & t<1.1)+ ...
%       cos(2*pi*50*t).*(t>= 0.2 & t< 1.4)+0.25*randn(size(t));
%   y = sin(2*pi*10*t).*(t>=0.6 & t<1.2)+...
%       sin(2*pi*50*t).*(t>= 0.4 & t<1.6)+ 0.35*randn(size(t));
%   wcoherence(x,y,1000)
%
%   % Example 2:
%   %   Plot the wavelet coherence between the El Nino time series and the
%   %   All Indian Average Rainfall Index. The data are sampled monthly.
%   %   Set the phase display threshold to 0.7. Specify the sampling
%   %   interval as 1/12 of a year to display the periods in years.
%   load ninoairdata;
%   wcoherence(nino,air,years(1/12),'phasedisplaythreshold',0.7);
%
%   See also cwtft, duration


narginchk(1,11);
nargoutchk(0,1);

%Check input vector size
nx = numel(x);
validateattributes(x,{'numeric'},{'real','finite'},'wcoherence','X');

% Form signals as row vectors
x = x(:)';


params = parseinputs(nx,varargin{:});

% Get number of voices per octave
nv = params.nv;

% If sampling frequency is specified, dt = 1/fs
if (isempty(params.fs) && isempty(params.Ts))
    % The default sampling interval is 1 for normalized frequency
    dt = params.dt;
    
elseif (~isempty(params.fs) && isempty(params.Ts))
    % Accept the sampling frequency in hertz
    fs = params.fs;
    dt = 1/fs;
elseif (isempty(params.fs) && ~isempty(params.Ts))
    % Get the dt and Units from the duration object
    [dt,~] = getDurationandUnits(params.Ts);
    
    
end

%Create scale vector for the CWT
s0 = 2*dt;
a0 = 2^(1/nv);
noct = params.numoct;
scales = s0*a0.^(0:noct*nv);
scales = scales(:);
invscales = 1./scales;
invscales = repmat(invscales,1,nx);

ns = params.numscalestosmooth;

wname = 'morl';
cwtx = cwtft({x,dt},'wavelet',wname,'scales',scales,'PadMode','symw');

cwtx.cfs = cwtx.cfs(:,1:nx);
cwtx.cfs1 = smoothCFS(invscales.*abs(cwtx.cfs).^2,scales,dt,ns);

end

function cfs = smoothCFS(cfs,scales,dt,ns)
N = size(cfs,2);
npad = 2.^nextpow2(N);
omega = 1:fix(npad/2);
omega = omega.*((2*pi)/npad);
omega = [0., omega, -omega(fix((npad-1)/2):-1:1)];

% Normalize scales by DT because we are not including DT in the
% angular frequencies here. The smoothing is done by multiplication in
% the Fourier domain
normscales = scales./dt;
for kk = 1:size(cfs,1)
    F = exp(-0.25*(normscales(kk)^2)*omega.^2);
    smooth = ifft(F.*fft(cfs(kk,:),npad));
    cfs(kk,:)=smooth(1:N);
end
% Convolve the coefficients with a moving average smoothing filter across
% scales
H = 1/ns*ones(ns,1);
cfs = conv2(cfs,H,'same');
end
%------------------------------------------------------------------------
function params = parseinputs(N,varargin)
% Set up defaults
params.fs = [];
params.dt = 1;
params.Ts = [];
params.sampinterval = false;
params.engunitflag = true;
params.normalizedfreq = true;
params.nv = 12;
params.numscalestosmooth = 12;
maxnumoctaves = floor(log2(N))-1;
params.numoct = maxnumoctaves;
params.wav = 'morl';

% Error out if there are any calendar duration objects
tfcalendarDuration = cellfun(@iscalendarduration,varargin);
if any(tfcalendarDuration)
    error(message('Wavelet:FunctionInput:CalendarDurationSupport'));
end

tfsampinterval = cellfun(@isduration,varargin);

if (any(tfsampinterval) && nnz(tfsampinterval) == 1)
    params.sampinterval = true;
    params.Ts = varargin{tfsampinterval>0};
    if (numel(params.Ts) ~= 1 ) || params.Ts <= 0 || isempty(params.Ts)
        error(message('Wavelet:FunctionInput:PositiveScalarDuration'));
    end
    
    params.engunitflag = false;
    params.normalizedfreq = false;
    varargin(tfsampinterval) = [];
end

params.mincoherence = 0.5;
tfvoices = find(strncmpi(varargin,'voicesperoctave',1));
if any(tfvoices)
    
    params.nv = varargin{tfvoices+1};
    validateattributes(params.nv,{'numeric'},{'positive','integer',...
        'scalar','>=',10,'<=',32},'wcoherence','VoicesPerOctave');
    varargin(tfvoices:tfvoices+1) = [];
end

tfnumoctaves = find(strncmpi(varargin,'numoctaves',4));

if any(tfnumoctaves)
    params.numoct = varargin{tfnumoctaves+1};
    validateattributes(params.numoct,{'numeric'},{'positive','integer',...
        '<=',maxnumoctaves},'wcoherence','NumOctaves');
    varargin(tfnumoctaves:tfnumoctaves+1) = [];
end

%The number of scales to smooth defaults to the number of voices per
%octave
params.numscalestosmooth = params.nv;

tfmincoherence = find(strncmpi(varargin,'phasedisplaythreshold',1));

if any(tfmincoherence)
    params.mincoherence = varargin{tfmincoherence+1};
    validateattributes(params.mincoherence,{'numeric'},{'scalar','>=',0,...
        '<=',1},'wcoherence','PhaseDisplayThreshold');
    varargin(tfmincoherence:tfmincoherence+1) = [];
end

maxsmooth = floor((params.nv*params.numoct+1)/2);
tfnumscalestosmooth = find(strncmpi(varargin,'numscalestosmooth',4));

if any(tfnumscalestosmooth)
    params.numscalestosmooth = varargin{tfnumscalestosmooth+1};
    validateattributes(params.numscalestosmooth,{'numeric'},{'positive',...
        'integer','scalar','<=',maxsmooth},'wcoherence','NumScalesToSmooth');
    varargin(tfnumscalestosmooth:tfnumscalestosmooth+1) = [];
end

% Only scalar left must be sampling frequency
tfsampfreq = cellfun(@(x) (isscalar(x) && isnumeric(x)),varargin);

if (any(tfsampfreq) && (nnz(tfsampfreq) == 1) && ~params.sampinterval)
    params.fs = varargin{tfsampfreq};
    validateattributes(params.fs,{'numeric'},{'positive'},'wcoherence','Fs');
    params.normalizedfreq = false;
    params.engunits = true;
elseif any(tfsampfreq) && params.sampinterval
    error(message('Wavelet:FunctionInput:SamplingIntervalOrDuration'));
elseif nnz(tfsampfreq)>1
    error(message('Wavelet:FunctionInput:Invalid_ScalNum'));
end
end
