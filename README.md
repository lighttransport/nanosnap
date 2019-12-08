# NanoSNAP, Nanoscale Signal, Noise and Audio Processing library in C++11

NanoSNAP is a small and portable signal, audio and noise processing library in C++11. 🤞
NanoSNAP depends only on C++11 STL.

## Usage

* For running TTS(Text-to-speech) and ASR(Automatic Speech Recognition) on C++ Embedded device.
* Image processing with neural netowork inference on C++ and Embedded device.
* Implement audio and speech feature(e.g. using `rfft`, `mfcc` `stft`, `istft`) on your C++ machine learning library.

## Install and integration

Simply copy `include` and `src` folder to your platform.

## Requirements

* CMake(for building examples and tests, and build NanoSNAP as submodules)
* C++11 compiler

## Supported platform

* [x] Windows
* [x] Linux
* [x] macOS
* [ ] Android(not tested yet, but should work)
* [ ] Raspberry Pi(not tested yet, but should work)
* [ ] RISC-V(not tested yet, but should work)

## Build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

### Build on Visual Studio 2017

```
> vcsetup.bat
```

Open `build/nanosnap.sln` and build it.

### Build and running tests

```
$ mkdir build
$ cd build
$ cmake -DNANOSNAP_ENABLE_TESTS=On ..
$ make
$ ./bin/test_nanosnap
```

### Compiler macros

* `NANOSNAP_NO_STDIO` Disable IO. e.g. `wav_read` is not available. This feature is useful when you want to use NanoSNAP in Android or embedded devices.


## API design

NanoSNAP takes raw pointer for input array values followin its length information(or shape information).

```c++
bool proc(const float *input, int n);
```

Output array is usually `std::vector` type so that NanoSNAP can allocate buffer for output internally.
Output array is a functiona argument when a function needs to return the status.

```c++
bool proc(int n, std::vector<float> *output);
```

Otherwise, output array is a return value.

```c++
std::vector<float> proc(int n);
```

### Internal state.

All API does not contain its internal state.

### Multithreading

NanoSNAP API does not ensure MT-safe.

### CMake option for developers

* `-DSANITIZE_ADDRESS=On` : Enable ASan.

## Data layout of array

NanoSNAP process 2D and higher ND array data as 1D flattened array.

The ordering of array data follows C language(This is same behavior in `numpy` array in C mode). For example, `img[H][W]` has `W` pixels in width(colums) , `H` pixels in height(rows).

```
-> memory address increases

+-----------+-----------+     +-------------+-----------+     +---------------+
| img[0][0] | img[0][1] | ... | img[0][W-1] | img[1][0] | ... | img[H-1][W-1] |
+-----------+-----------+     +-------------+-----------+     +---------------+
```

In contrary to `numpy` or vision/ML library such like OpenCV, The notation of dimensional arguments for a function signature starts from inner most dimension(right-most array dim). This is rather common notation in C language and graphics community. i.e,

```
// `output` has the shape of [h][w]
void create_image(size_t w, size_t h, float *output);

// `output` has the shape of [d][h][w]
void create_3d_tensor(size_t w, size_t h, size_t d, float *output);

// `input` has the shape of [nrows][nframes].
void rfft(const float *inout, size_t nframes, size_t nrows, ...);
```

## Features

### Numpy

| NanoSNAP               | Description                                        | Python equivalent                    |
| ---------------------- | -------------------------------------------------- | ------------------------------------ |
| `reshape_with_strides` | Create an array with the given shape and strides.  | `numpy.lib.stride_tricks.as_strided` |
| `convolve`             | 1D convolution                                     | `numpy.convolve`                     |
| `loadtxt`              | Load 1D or 2D array                                | `numpy.loadtxt`                      |
| `savetxt`              | Save 1D or 2D array                                | `numpy.savetxt`                      |

### Random number generation

| NanoSNAP               | Description            | Python equivalent                  |
| ---------------------- | ---------------------- | ---------------------------------- |
| `random_uniform`       | Uniform random number  | `numpy.random.rand`                |
| `random_shuffle`       | Randomly shuffle array | `numpy.random.shuffle`             |


### FFT

| NanoSNAP               | Description                  | Python equivalent                   |
| ---------------------- | ---------------------------- | ----------------------------------- |
| `rfft`                 | Real 1D FFT                  | `numpy.fft.rfft`                    |
| `ifft`                 | Inverse Complex FFT          | `numpy.fft.ifft`                    |

### Scipy

| NanoSNAP               | Description                                                 | Python equivalent                   |
| ---------------------- | ----------------------------------------------------------- | ----------------------------------- |
| `lfilter`              | Filter data along one-dimension with an IIR or FIR filter.  | `scipy.signal.lfilt`                |
| `medfilt`              | Median filter                                               | `scipy.signal.medfilt`              |
| `wav_read`             | Read .WAV file                                              | `scipy.io.wavfile.read`             |
| `wav_write`            | Save .WAV file                                              | `scipy.io.wavfile.write`            |

### Python speech features

| NanoSNAP               | Description                                       | Python equivalent                   |
| ---------------------- | ------------------------------------------------- | ----------------------------------- |
| `mel2hz`               | Mel to Hz                                         | `mel2hz`                            |
| `hz2mel`               | Hz to Mel                                         | `hz2mel`                            |
| `lifter`               | Apply a cepstral lifter the the matrix of cepstra | `lifter`                            |

#### Work in progress

| NanoSNAP               | Description                                       | Python equivalent                   |
| ---------------------- | ------------------------------------------------- | ----------------------------------- |
| `mfcc`                 | Mel Frequency Cepstral Coefficients               | `mfcc`                              |
| `fbank`                | Filterbank Energies                               | `fbank`                             |
| `logfbank`             | Log Filterbank Energies                           | `logfbank`                          |
| `ssc`                  | Spectral Subband Centroids                        | `ssc`                              |


### Librosa

| NanoSNAP               | Description                                                            | Python equivalent                   |
| ---------------------- | ---------------------------------------------------------------------- | ----------------------------------- |
| `stft`                 | Short Term Fourier Transform                                           | `librosa.stft`                      |
| `istft`                | Inverse STFT                                                           | `librosa.istft`                     |
| `mel`                  | Create a Filterbank matrix to combine FFT bins into Mel-frequency bins | `librosa.filters.mel`                     |

### Image

| NanoSNAP               | Description                  | Python equivalent                   |
| ---------------------- | ---------------------------- | ----------------------------------- |
| `resize_bilinear`      | Resize image with bilinear   | `cv2.resize_image`                  |
| `imread`               | Load LDR image               | `cv2.imread`                        |
| `imsave`               | Save image as LDR format     | `cv2.imsave`                        |

## limited support

* get_window : 'hann' only. `scipy.signal.get_window` equivalent.

## TODO

* [ ] Image reszier.
* [ ] Better error handling(report error message)
* [ ] Multithreading using C++11 `thread`.
  * [ ] Use `StackVector` as much as possible.
* [ ] Read/write WAV from/to buffer(memory)
* [ ] Integrate with `NanoNumCp`
* FFT
  * [ ] Implement more FFT functions defined in `scipy.fft`.
  * [ ] 2D FFT
  * [ ] Replace pocketfft with muFFT https://github.com/Themaister/muFFT or our own C++11 FFT routuine.
* [ ] Port more functions in `python_speech_features`
* [ ] Implement more speech features implemented in sox, librosa, etc.
* [ ] Plot and save figure/image in JPG/PNG/EXR

## Developer note

### Adding tests

* Write testvector generator and put it to `tests/gen/`
  * Generate testvector file(`.inc`)
* Add .cc to `tests`. Add it to CMakeLists.txt.

## License

NanoSNAP is licensed under MIT license.

### Third party licenses.

* stack_vector.h : Copyright (c) 2006-2008 The Chromium Authors. All rights reserved. Use of this source code is governed by a BSD-style license.
* doctest : The MIT License (MIT). Copyright (c) 2016-2019 Viktor Kirilov
* dr_wav : Public domain or MIT-0. By David Reid.
* python_speech_features : The MIT License (MIT). Copyright (c) 2013 James Lyons. https://github.com/jameslyons/python_speech_features
* pocketfft : FFT library used in numpy. Copyright (C) 2004-2018 Max-Planck-Society. 3-clause BSD-tyle license. https://gitlab.mpcdf.mpg.de/mtr/pocketfft
* c_speech_features : Copyright (c) 2017 Chris Lord. MIT license. https://github.com/Cwiiis/c_speech_features
* STB image : Public domain. https://github.com/nothings/stb
* sRGB transform : Copyright (c) 2017 Project Nayuki. (MIT License) https://www.nayuki.io/page/srgb-transform-library
