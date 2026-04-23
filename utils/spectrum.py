import numpy as np
import scipy.signal as signal


def cepstrum(sig, fs=1.0):
    freq = np.fft.rfftfreq(sig.shape[0], 1 / fs)
    coef = np.fft.rfft(sig)

    log_abs_X = np.log(np.abs(coef))
    log_abs_X -= np.mean(
        log_abs_X)  # Remove DC component for cepstrum calculation

    cepstrum = np.fft.rfft(log_abs_X)
    df = freq[1] - freq[0]
    quefrency = np.fft.rfftfreq(log_abs_X.size, df)

    return cepstrum, quefrency
