# NanoSNP, Nanoscale Signal, Audio and Noise Processing library in C++11

NanoSNP is a small and portable signal, audio and noise processing library in C++11.

## Install

Simply copy `include` and `src` folder to your platform.

## Requirements

* CMake(for building examples and tests, and build NanoSNP as submodules)
* C++11 compiler

## Supported platform

* [x] Windows
* [x] Linux
* [x] macOS
* [ ] Android
* [ ] Raspberry Pi

## Build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Features

| NanoSNP                | Python equivalent                   |
| ---------------------- | ----------------------------------- |
| `medfilt`              | `scipy.signal.medfilt`              |
| `wav_read`             | `scipy.io.wavfile.read`             |
| `wav_read_from_buffer` | Read WAV from buffer(memory)        |
| `wav_write`            | `scipy.io.wavfile.write`            |
| `wav_write_to_buffer`  | Write WAV to buffer(memory)         |
|                        | `python_speech_features.fbank`      |
|                        | `python_speech_features.logfbank`   |
|                        | `python_speech_features.mfcc`       |
|                        | `python_speech_features.ssc`        |

## Compiler macros

* `NANOSNP_NO_STDIO` Disable IO. e.g. `wav_read` is not available. This feature is useful when you want to use NanoSNP in Android or embedded devices.

## TODO

* [ ] Write our own FFT routine.

## License

NanoSNP is licensed under MIT license.

### Third party licenses.

* doctest : The MIT License (MIT). Copyright (c) 2016-2019 Viktor Kirilov
* dr_wav : Public domain or MIT-0. By David Reid.
* fft2d : Very permissive license. See `src/fft2d/readme.txt` for details. Copyright(C) 1996-2001 Takuya OOURA

* T.B.W.
