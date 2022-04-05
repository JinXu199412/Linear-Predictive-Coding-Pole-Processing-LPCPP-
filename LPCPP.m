%% Linear Predictive Coding Pole Processing
% Linear Predictive Coding Pole Processing (LPCPP) is parameterised
% time-frequency method which can be used to identify and track the
% dominant frequency change in real-time (i.e. the term 'real-time' refers
% to every sampling instant of a discrete signal).

% The LPCPP method further processes LPC (Linear Predictive Coding) poles
% after the LPC method to produce a series of reduced-order filter
% transform functions to realise the dominant frequency estimation.
% Specifically, there are three steps:
% 1. Categorise LPC poles into dominant poles and non-dominant poles.
% 2. Identify the associated poles of each dominant pole from non-dominant 
% poles.
% 3. The dominant poles and their corresponding associated non-dominant 
% pole(s) are used to form a series of reduced-order transform functions.
% Details of the LPCPP method can refer to our publications.

% Author(s): Jin Xu, Mark Davis, Ruairi de Frein, 2022-03-25
% Publication(s): 
% [1] Xu, J., Davis, M. and de Frein, R., 2021. An LPC pole processing
% method for enhancing the identification of dominant spectral features. 
% Electronics Letters, 57(18), pp.708-710. 
% (https://doi.org/10.1049/ell2.12226)
% [2] Xu, J., Davis, M. and de Frein, R., 2022. Dominant frequency
% component tracking of noisy time-varying signals using the linear 
% predictive coding pole processing method. Electronics Letters, 58(2),
% pp.79-81. (https://doi.org/10.1049/ell2.12362)
% Acknowledgements: 
% This work is from a PhD thesis at Technological University Dublin.

%% LPCPP MATLAB Function
function [fre_estimates] = LPCPP(signal,fs,P,beta,lambda)
% fre_estimates = LPCPP(signal,P,fs,beta,lambda) returns the dominant
% frequency estimates of a signal. signal is a vector. P is the number of
% filter orders. fs is the signal sampling frequency. beta is a threshold
% value for identifying the dominant pole(s) and its value range is (0,1).
% lambda is a threshold value for identify the assoicated pole(s) for each
% dominant pole and the choice of the lambda value is recommended to be 
% between [5,10].

%---------------------------------
% Check inputs/outputs
if nargin < 1 
    msg = 'Function input cannot be empty!';
    error(msg)
end
[m,n] = size(signal);
if (m>1) && (n>1)
    msg = 'Input signal must be a vector!';
    error(msg)
end
if (m>1) && (n==1)
	signal = signal';
	[~,n] = size(signal);
end
if nargin < 2 || isempty(fs)
    msg = 'fs input cannot be empty!';
    error(msg)
end   
if nargin < 3 || isempty(P)
    P=10;
end
if nargin < 4 || isempty(beta)
    beta = 0.5;
end
if nargin < 5 || isempty(lambda)
    lambda = 10;
end

%---------------------------------
% Compute the LPC coefficients

% LPCPP frequency estimate precision is 0.01
lpc_f = 0:0.01:fs/2;
f_L = 2*length(lpc_f);
hf_l = floor(f_L/2);
% Compute autocorrelation vector or matrix
X = fft(signal,f_L);
R = ifft(abs(X).^2);
R = R(1:P+1);
acf = R';
% Build Toeplitz matrix
RM = toeplitz(acf(1:P), conj(acf(1:P)));
% Obtain LPC coefficients
LPC_a = RM \ -acf(2:end);
LPC_a = [1, LPC_a'];

%---------------------------------
% Categorise LPC poles into dominant pole and non-dominant poles

% Calculate the LPC roots (poles)
all_rts = roots(LPC_a);
% Only consider the poles with non-negative imaginary parts
rts = all_rts(imag(all_rts)>0);
% Magnitude of the LPC poles
ampt_poles = abs(rts);
[ampt_poles,I] = sort(ampt_poles,'descend');
rts = rts(I);
angz = atan2(imag(rts),real(rts));
poles = angz.*(fs/(2*pi));
% Single-pole transfer function for each LPC pole
all_Ho = zeros(length(rts),hf_l);
for d_i = 1:length(rts)
    Ho =fft([1 -rts(d_i)],f_L);
    Ho = 1./Ho;
    Ho = Ho(1:hf_l);
    Ho = abs(Ho).^2;
    all_Ho(d_i,:)=Ho;
end
% Enhancement function for LPC poles magnitude
enhanced_m = 1./(1-ampt_poles);
% Identify the dominant poles
sumall_em = sum(enhanced_m);
i = 1;
manitude_sum =0;
while manitude_sum < beta*sumall_em
    manitude_sum = enhanced_m(i)+manitude_sum;
    i = i+1;
end
dominant_l = i-1;
dominant_poles = poles(1:dominant_l,1)';
dominant_spectrum = all_Ho(1:dominant_l,:);
nondominant_pole = poles(dominant_l+1:end,1)';
nondominant_spectrum=all_Ho(dominant_l+1:end,:);

%---------------------------------
% Identify the associated poles of each dominant pole

ass_l = length(nondominant_pole);
ass_domin_distance = zeros(i-1,ass_l);
fre_estimates = zeros(1,dominant_l);
for d_i = 1:dominant_l
    % Identify the associated poles
    ass_domin_distance(d_i,:)=abs(nondominant_pole-dominant_poles(d_i));
    d_a_index = ass_domin_distance(d_i,:) < lambda;
    if sum(d_a_index)==1
        ghost_H = nondominant_spectrum(d_a_index,:);
    else
        ghost_H = prod(nondominant_spectrum(d_a_index,:));
    end
 
%---------------------------------
% Form a series of reduced-order transform functions and give the results

    % Form a reduced-order transform function
    ghost_H = prod([ghost_H;dominant_spectrum(d_i,:)]);
    [~,I]=max(ghost_H);
    % Obtain frequency estimate
    peak_f = lpc_f(I);
    % Set outputs
    fre_estimates(1,d_i)= peak_f;
end

end

