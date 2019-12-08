python gen_medfilt1.py > ../testvector/medfilt1.inc
python gen_lifter.py > ../testvector/lifter.inc
python gen_rfft.py > ../testvector/rfft.inc
python gen_ifft.py > ../testvector/ifft.inc
python gen_librosa_stft.py > ../testvector/librosa_stft.inc
python gen_librosa_istft.py > ../testvector/librosa_istft.inc
python gen_librosa_filters_mel.py > ../testvector/librosa_filters_mel.inc
python gen_random_uniform.py > ../testvector/random_uniform.inc
python gen_random_normal.py > ../testvector/random_normal.inc
python gen_convolve_full.py > ../testvector/convolve_full.inc
python gen_convolve_same.py > ../testvector/convolve_same.inc
python gen_convolve_valid.py > ../testvector/convolve_valid.inc

python gen_signal_get_window_hann.py > ../testvector/signal_get_window_hann.inc
