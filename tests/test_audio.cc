#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <iostream>

TEST_CASE("wav_read_16bit") {
  const std::string &filename = "../tests/assets/8k16bitpcm.wav";

  uint32_t rate;
  std::string dtype;
  uint32_t channels;
  uint64_t samples;
  std::vector<uint8_t> data;

  std::string err;

  bool ret = nanosnap::wav_read(filename, &rate, &dtype, &channels, &samples, &data, &err);

  if (!err.empty())
  {
    std::cerr << "ERR: " << err << std::endl;
  }

  std::cout << "rate = " << rate << ", dtype = " << dtype << "\n";
  std::cout << "channels = " << channels << ", samples = " << samples << "\n";
  std::cout << "data.size = " << data.size() << "\n";


  CHECK(ret == true);
  CHECK(err.empty());
  CHECK(channels == 1);
  CHECK(dtype.compare("int16") == 0);
}

TEST_CASE("wav_read_8bit") {
  const std::string &filename = "../tests/assets/8k8bitpcm.wav";

  uint32_t rate;
  std::string dtype;
  uint32_t channels;
  uint64_t samples;
  std::vector<uint8_t> data;

  std::string err;

  bool ret = nanosnap::wav_read(filename, &rate, &dtype, &channels, &samples, &data, &err);

  if (!err.empty())
  {
    std::cerr << "ERR: " << err << std::endl;
  }

  std::cout << "rate = " << rate << ", dtype = " << dtype << "\n";
  std::cout << "channels = " << channels << ", samples = " << samples << "\n";
  std::cout << "data.size = " << data.size() << "\n";


  CHECK(ret == true);
  CHECK(err.empty());
  CHECK(channels == 1);
  CHECK(dtype.compare("uint8") == 0);
}
