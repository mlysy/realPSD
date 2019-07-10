#' @title Generate_FreqNoise_FromAPSD
#' @descrb Generates a time-domain noise array of total duration T (s) with a sampling frequency SF (Hz). The corresponding frequency-domain noise is generated assuming independent Gaussian distributed amplitudes for both real and imaginary components of the Fourier basis random variable pairs. The iFFT function is used to generate the corresponding time-domain array. The frequency span of APSD should exceed the span of the desired noise signal bandwidth. Otherwise, an error occurs during interpolation.
#' 
#' @param xAPSD Frequency vector of yAPSD (Hz) : [1 2 3 ... N/2]*FR this vector need NOT be linearly spaced.
#' @param yAPSD Average power spectral density (u^2/Hz).
#' @param SF    Sampling frequency (Hz).
#' @param T     Total duration (s).
#' @param Q     Quality factor to scale sine-wave amplitude by.
#' @param displ 1: display information, 0, don't display information.
#' @return A list of xTime = time (s) and yTime = amplitude (u)
generate_freqNoise_fromAPSD <- function(xAPSD, yAPSD, 
  SF, T, Q, CORRUPT, displ, varargin) {
  yPSD_noise = NA
  T = 2*T # doubling the time, to throw out half the data afterwards 
  xAPSD = c(xAPSD[1]/2, xAPSD) # extending the frequency range down to 1/T/2
  yAPSD = c(yAPSD[1], yAPSD)
  # --- Override resonance frequency if needed------------------------------
  if(length(varargin) < 1)
    f0_s = 33553
  else
    # f0_s = varargin{1}
    f0_s = varargin
  # ERROR PREVENTION - Extend the frequency range if necessary
  if(xAPSD[1] > 1/T) {
    print('!!! Cannot generate time signal longer than the 1/FR, where FR = the frequency resolution of the PSD.')
    print('      - Extending white noise to low frequency to correct the problem.');
    xAPSD = c(1/T, xAPSD)
    yAPSD = c(yAPSD[1], yAPSD)
  }

  if(xAPSD[length(xAPSD)]*2 < SF) {
    print('!!! Cannot generate time signal with Nyquist frequency higher than the maximum PSD frequency.');
    print('      - Extending white noise to high frequency to correct the problem.');
    xAPSD = c(xAPSD, SF)
    yAPSD = c(yAPSD, yAPSD[length(yAPSD)])
  }

  # # Determine total size N of the output yTime vector
  # N = floor(T*SF/2)*2   # round to nearest even integer
  # T = N/SF              # redefine the total duration after rounding

  # # Construct the desired xPSD vector, corresponding to xTime
  # FR = 1/T              # frequency resolution
  # xPSD = [FR:FR:SF/2]   # where SF/2 is the Nyquist frequency

#   # -------------------- ongoing -----------------
#   if(!is.null('displ', 'var') && displ==1){
#         print('')
#         print(['Generating time signal from APSD. N=' num2str(N/2) ' ; T=' num2str(T/2) 's'])
#   }
#   # Interpolate yAPSD onto the smaller frequency basis xPSD
#   yAPSD_inter = interp1(xAPSD, yAPSD, xPSD);

#   snan = sum(isnan(yAPSD_inter));
#   if snan>1
#       disp(['WARNING: NAN values during interpolation... '])
#   end

# % Sample the yASD_inter to obtain generated noise signal in the frequency-domain: yFFT
# % Multiply each yAPSD-inter value by a complex Gaussian random variable, with components ~N(0,1)
# K = length(yAPSD_inter);
# a = randn(1,K);
# b = randn(1,K);
# yFFT = sqrt(yAPSD_inter/2).*(a+b*1i);

# yFFT(1:end-1) = yFFT(1:end-1)/sqrt(2); % divide by sqrt(2) because the components (except Nyquist) will be unfolded into negative frequencies.
# xPSD = linspace(FR, SF/2, N/2);

# % Add sine-wave PSD if flag is true
# if CORRUPT == 1
#     f_noise = f0_s + 10*randn; % (Hz), with random jitter
#     A_noise = (Q^0.5).*3.5e3; % amplitude
#     [~, ~, xFreq, yFreq] = noise_sin(f_noise, A_noise, T, SF, 0);
#     assert(sum(xPSD ~= xFreq) == 0)
#     yFreq = yFreq./(sqrt(FR)*N); % scale to match yFFT scaling
#     yFFT_noise = yFFT + yFreq; % adding together (both unfolded)
# else
#     yFFT_noise = NaN;
# end

# % ----No noise added------
# yFFT(end) = real(yFFT(end)); % Nyquist freq should be real for real signal
# yPSD = abs(yFFT).^2; % Take the squared magnitude at each frequency
# % Make the PSD single-sided by doubling everything except the Nyquist frequency.
# % This corrects for the loss in amplitude caused by dropping negative frequencies
# % the Nyquist frequency at N/2 doesn't undergo folding, and therefore is NOT doubled.
# yPSD(1:N/2-1) = 2.*yPSD(1:N/2-1);

# % ----Sine-noise added------
# if CORRUPT == 1
#     yFFT_noise(end) = real(yFFT_noise(end)); % Nyquist freq should be real for real signal
#     yPSD_noise = abs(yFFT_noise).^2; % Take the squared magnitude at each frequency
#     % Make the PSD single-sided by doubling everything except the Nyquist frequency.
#     % This corrects for the loss in amplitude caused by dropping negative frequencies
#     % the Nyquist frequency at N/2 doesn't undergo folding, and therefore is NOT doubled.
#     yPSD_noise(1:N/2-1) = 2.*yPSD_noise(1:N/2-1);
#     if all(size(yPSD_noise) ~= size(xPSD)), xPSD = xPSD'; end
# end

# if all(size(yPSD) ~= size(xPSD)), xPSD = xPSD'; end

# if displ==1
#     disp(['Completed in ' num2str(t) ' seconds'])
# end
}

