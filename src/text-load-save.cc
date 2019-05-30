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
