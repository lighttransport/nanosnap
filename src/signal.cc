#include "nanosnap/nanosnap.h"

#include "stack_vector.h"

#include <cmath>
#include <cstring>
#include <memory>

//#include <iostream> // dbg

namespace nanosnap {

namespace {

// Find median value.
// `a` will be modified.
static inline float find_median(const size_t n, float *a) {
  // Sort `a`.
  // Simple O(n^2) selection sort.
  //
  for (size_t j = 0; j < n - 1; j++) {
    size_t iMin = j;

    /* test against elements after j to find the smallest */

    for (size_t i = j + 1; i < n; i++)

    {
      /* if this element is less, then it is the new minimum */

      if (a[i] < a[iMin])

      {
        /* found new minimum; remember its index */

        iMin = i;
      }
    }

    if (iMin != j) {
      std::swap(a[j], a[iMin]);
    }
  }

  // We know the window size is odd, thus a[n / 2] is the median value.
  return a[(n >> 1)];
}

}  // namespace

bool medfilt1(const float *x, const size_t n, const int k, std::vector<float> *y,
              bool include_nan, bool padding) {
  if ((k % 2) != 1) {
    // window size must be odd.
    return false;
  }

  // TODO(LTE): Handle these parameters.
  (void)padding;

  StackVector<float, 64> buf;
  buf->resize(size_t(k * k));

  y->resize(n);

  // TODO(LTE): Test if k < 2

  int m = (k - 1) / 2;
  for (ssize_t i = 0; i < ssize_t(n); i++) {
    buf->clear();
    for (ssize_t r = -m; r <= m; r++) {
      ssize_t j = i + r;

      float v;
      bool out_of_bounds = (j < 0) || (j >= ssize_t(n));
      if (out_of_bounds) {
        v = 0.0f;
      } else {
        v = x[size_t(j)];
      }

      if (!include_nan) {
        if (std::isnan(v)) {
          v = 0.0f;
        }
      }

      buf->push_back(v);
    }

    (*y)[size_t(i)] = find_median(size_t(k), buf->data());
  }

  return true;
}


bool convolve(const float *_a, const size_t _n, const float *_v,
              const size_t _m, std::vector<float> *output, const int mode) {
  size_t m = _m, n = _n;
  const float *a = _a;
  const float *v = _v;

  if (_m > _n) {
    // swap
    a = _v;
    n = _m;

    v = _a;
    m = _n;
  }

  if ((n == 0) || (m == 0)) {
    return false;
  }

  if (mode == 0) {
    // mode 'full'
    //
    // len(a) = 5
    // len(v) = 3
    // output = 5 + 3 - 1 = 7
    //
    //                 +----+----+----+----+----+
    //                 | a0 | a1 | a2 | a3 | a4 |
    //                 +----+----+----+----+----+
    //
    //             ....+----+
    //  n0     v2 | v1 | v0 |
    //             ....+----+
    //
    //             ....+----+----+
    //  n1         v2  | v1 | v0 |
    //             ....+----+----+
    //
    //                 +----+----+----+
    //  n2             | v2 | v1 | v0 |
    //                 +----+----+----+
    //
    //                                     +----+....
    //  n6                                 | v2 | v1 | v0
    //                                     +----+....
    //

    size_t out_length = n + m + 1;
    output->resize(out_length);

    // TODO(LTE): optimize loop
    for (int n_idx = 0; n_idx < int(out_length); n_idx++) {
      float sum = 0.0f;
      for (int m_idx = -(int(m) - 1); m_idx <= 0; m_idx++) {
        int a_idx = n_idx + m_idx;
        int v_idx = -m_idx;
        // std::cout << "n = " << n_idx << ", m = " << m_idx << "a = " << a_idx
        // << ", v = " << v_idx << "\n";
        if ((a_idx < 0) || (a_idx >= int(n))) {
          continue;
        }
        if ((v_idx < 0) || (v_idx >= int(m))) {
          continue;
        }
        // std::cout << "got it\n";
        sum += a[size_t(a_idx)] * v[size_t(v_idx)];
      }
      // std::cout << "sum[" << n_idx << "] = " << sum << std::endl;
      (*output)[size_t(n_idx)] = sum;
    }
  } else if (mode == 1) {
    // mode 'same'
    //
    // len(a) = 5
    // len(v) = 3
    // output = 5 (max(M, N))
    //
    //                 +----+----+----+----+----+
    //                 | a0 | a1 | a2 | a3 | a4 |
    //                 +----+----+----+----+----+
    //
    //             ....+----+----+
    //  n0         v2  | v1 | v0 |
    //             ....+----+----+
    //
    //                 +----+----+----+
    //  n2             | v2 | v1 | v0 |
    //                 +----+----+----+
    //
    //                                +----+----+....
    //  n4                            | v2 | v1 | v0
    //                                +----+----+....
    //

    size_t out_length = n;
    output->resize(out_length);

    // TODO(LTE): Verify this offset calculation is correct.
    int a_offset = int(m) / 2;

    // TODO(LTE): optimize loop
    for (int n_idx = 0; n_idx < int(out_length); n_idx++) {
      float sum = 0.0f;
      for (int m_idx = -(int(m) - 1); m_idx <= 0; m_idx++) {
        int a_idx = n_idx + m_idx + a_offset;
        int v_idx = -m_idx;
        // std::cout << "n = " << n_idx << ", m = " << m_idx << "a = " << a_idx
        // << ", v = " << v_idx << "\n";
        if ((a_idx < 0) || (a_idx >= int(n))) {
          continue;
        }
        if ((v_idx < 0) || (v_idx >= int(m))) {
          continue;
        }
        // std::cout << "got it\n";
        sum += a[size_t(a_idx)] * v[size_t(v_idx)];
      }
      // std::cout << "sum[" << n_idx << "] = " << sum << std::endl;
      (*output)[size_t(n_idx)] = sum;
    }

  } else if (mode == 2) {
    // mode 'valid'
    //
    // len(a) = 5
    // len(v) = 3
    // output = 5 - 3 + 1
    //
    //                 +----+----+----+----+----+
    //                 | a0 | a1 | a2 | a3 | a4 |
    //                 +----+----+----+----+----+
    //
    //                 +----+----+----+
    //  n0             | v2 | v1 | v0 |
    //                 +----+----+----+
    //
    //                      +----+----+----+
    //  n1                  | v2 | v1 | v0 |
    //                      +----+----+----+
    //
    //                           +----+----+----+
    //  n2                       | v2 | v1 | v0 |
    //                           +----+----+----+
    //
    //

    size_t out_length = n - m + 1;
    output->resize(out_length);

    int a_offset = int(m) - 1;

    // TODO(LTE): optimize loop
    for (int n_idx = 0; n_idx < int(out_length); n_idx++) {
      float sum = 0.0f;
      for (int m_idx = -(int(m) - 1); m_idx <= 0; m_idx++) {
        int a_idx = n_idx + m_idx + a_offset;
        int v_idx = -m_idx;
        // std::cout << "n = " << n_idx << ", m = " << m_idx << "a = " << a_idx
        // << ", v = " << v_idx << "\n";
        if ((a_idx < 0) || (a_idx >= int(n))) {
          continue;
        }
        if ((v_idx < 0) || (v_idx >= int(m))) {
          continue;
        }
        // std::cout << "got it\n";
        sum += a[size_t(a_idx)] * v[size_t(v_idx)];
      }
      // std::cout << "sum[" << n_idx << "] = " << sum << std::endl;
      (*output)[size_t(n_idx)] = sum;
    }

  } else {
    return false;
  }

  return true;
}

bool lfilter(const float *b, const size_t nb, const float *a, const size_t na,
  const float *x, const size_t nx, const size_t mx)
{
  std::vector<float> b_normalized(nb);
  memcpy(b_normalized.data(), b, sizeof(float) * nb);

  // normalize `b' by a[0]
  for (size_t i = 0; i < nb; i++) {
    b_normalized[i] /= a[0];
  }

  // TODO(LTE): Implement
  (void)na;

/*
        ind = out_full.ndim * [slice(None)]
        if zi is not None:
            ind[axis] = slice(zi.shape[axis])
            out_full[ind] += zi

        ind[axis] = slice(out_full.shape[axis] - len(b) + 1)
out = out_full[ind]
*/

  // out_full = np.apply_along_axis(lambda y: np.convolve(b, y), axis, x)
  // Apply convolve for each row(axis = -1 behavior in scipy.signal.lfilter).
  for (size_t j = 0; j < mx; j++) {
    std::vector<float> output;
    bool ret = convolve(b, nb, &x[j * nx], nx, &output, /* mode */0);
    if (!ret) {
      return false;
    }
  }

  return false;

}

}  // namespace nanosna
