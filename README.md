# NanoSNAP, Nanoscale Signal, Noise and Audio Processing library in C++11

NanoSNAP is a small and portable signal, audio and noise processing library in C++11. 🤞

## Install

Simply copy `include` and `src` folder to your platform.

## Requirements

* CMake(for building examples and tests, and build NanoSNAP as submodules)
* C++11 compiler

## Supported platform

* [x] Windows
* [x] Linux
* [x] macOS
* [ ] Android
* [ ] Raspberry Pi
* [ ] RISC-V

## Build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

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

In contrary to `numpy` or vision/ML community, The notation of dimensional arguments for a function signature starts from inner most dimension(right-most array dim). This is rather common notation in C language and graphics community. i.e,

```
// `output` has the shape of [h][w]
void create_image(size_t w, size_t h, float *output) {
}

// `output` has the shape of [d][h][w]
void create_3d_tensor(size_t w, size_t h, size_t d, float *output) {
}

// `input` has the shape of [nrows][nframes].
void rfft(size_t nframes, size_t nrows, const float *inout, ...) {
}
```

## Features

### Random number generation

| NanoSNAP               | Description            | Python equivalent                  |
| ---------------------- | ---------------------- | ---------------------------------- |
| `random_uniform`       | Uniform random number  | `numpy.random.rand`                |


### FFT

| NanoSNAP               | Description        | Python equivalent                   |
| ---------------------- | ------------------ | ----------------------------------- |
| `rfft`                 | Real 1D FFT        | `numpy.fft.rfft`                    |

### Scipy

| NanoSNAP               | Description        | Python equivalent                   |
| ---------------------- | ------------------ | ----------------------------------- |
| `medfilt`              | Median filter      | `scipy.signal.medfilt`              |
| `wav_read`             | Read .WAV file     | `scipy.io.wavfile.read`             |
| `wav_write`            | Save .WAV file     | `scipy.io.wavfile.write`            |

### Python speech features

| NanoSNAP               | Description        | Python equivalent                   |
| ---------------------- | ------------------ | ----------------------------------- |
| `mel2hz`               | Mel to Hz          | `mel2hz`                            |
| `hz2mel`               | Hz to Mel          | `hz2mel`                            |
| `lifter`               |                    | `lifter`                            |
| `fbank`                |                    | `fbank`                             |


## TODO

* [ ] Multithreading with C++11 `thread`.
  * [ ] Use `StackVector` as much as possible.
* [ ] Read/write WAV from buffer(memory)
* [ ] Integrate with `NanoNumCp`
* FFT
  * [ ] Implement more FFT functions defined in `scipy.fft`.
  * [ ] 2D FFT
* [ ] Port `python_speech_features`
  * [ ] `python_speech_features.fbank`
  * [ ] `python_speech_features.logfbank`
  * [ ] `python_speech_features.mfcc`
  * [ ] `python_speech_features.ssc`
* [ ] Write our own FFT routine.

## License

NanoSNAP is licensed under MIT license.

### Third party licenses.

* stack_vector.h : Copyright (c) 2006-2008 The Chromium Authors. All rights reserved. Use of this source code is governed by a BSD-style license.
* doctest : The MIT License (MIT). Copyright (c) 2016-2019 Viktor Kirilov
* dr_wav : Public domain or MIT-0. By David Reid.
* fft2d : Very permissive license. See `src/fft2d/readme.txt` for details. Copyright(C) 1996-2001 Takuya OOURA
* python_speech_features : The MIT License (MIT). Copyright (c) 2013 James Lyons. https://github.com/jameslyons/python_speech_features
* pocketfft : FFT library used in numpy. Copyright (C) 2004-2018 Max-Planck-Society. 3-clause BSD-tyle license. https://gitlab.mpcdf.mpg.de/mtr/pocketfft

