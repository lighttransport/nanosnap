/*
The MIT License (MIT)

Copyright (c) 2019 Light Transport Entertainment, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#include "nanosnap/nanosnap.h"

#define NANOCSV_IMPLEMENTATION
#include "nanocsv.h"

namespace nanosnap {

bool loadtxt(const std::string &filename, std::vector<float> *values, int *n, int *m, std::string *err)
{
  nanocsv::ParseOption option;
  std::string warn;

  nanocsv::CSV<float> csv;

  bool ret = nanocsv::ParseCSVFromFile(filename, option, &csv, &warn, err);

  if (!ret) {
    return false;
  }

  (*n) = int(csv.num_fields);
  (*m) = int(csv.num_records);
  (*values) = std::move(csv.values);

  return true;

}

bool savetxt(const std::string &filename, const float *values, const int n, const int m, std::string *err)
{
  std::ofstream ofs(filename);

  if (!ofs) {
    if (err) {
      (*err) += "Failed to open file " + filename + " to write.\n";
    }
    return false;
  }

  // use ' ' as delimiter
  char delimiter = ' ';
  // TODO(LTE): Set precision

  for (size_t j = 0; j < size_t(m); j++) {
    for (size_t i = 0; i < size_t(n); i++) {
      ofs << values[j * size_t(n) + i];
      if (i != (size_t(n) - 1)) {
        ofs << delimiter;
      }
    }
    ofs << "\n";
  }

  return true;

}


} // namespace nanosnap
