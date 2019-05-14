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

## Features

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

* [ ] Read/write WAV from buffer(memory)
* [ ] Integrate with NanoNumCp
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

