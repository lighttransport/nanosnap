#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"


#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "nanosnap/nanosnap.h"

#include <algorithm>

namespace nanosnap {

namespace {

std::string GetFileExtension(const std::string &filename) {
  if (filename.find_last_of(".") != std::string::npos)
    return filename.substr(filename.find_last_of(".") + 1);
  return "";
}

// chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
// Fast and efficient sRGB -> Linear conversion. Error is small(abs diff <= 0.0018)
inline float sRGBToLinearFast(const float &sRGB) {
  const float Linear = sRGB * (sRGB * (sRGB * 0.305306011f + 0.682171111f) + 0.012522878f);
  return Linear;
}

/*
 * sRGB transform (C++)
 *
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/srgb-transform-library
 */
const float SRGB_8BIT_TO_LINEAR_FLOAT[1 << 8] = {
	0.0f, 3.03527e-4f, 6.07054e-4f, 9.10581e-4f,
	0.001214108f, 0.001517635f, 0.001821162f, 0.0021246888f,
	0.002428216f, 0.002731743f, 0.00303527f, 0.0033465358f,
	0.0036765074f, 0.004024717f, 0.004391442f, 0.0047769537f,
	0.005181517f, 0.005605392f, 0.0060488335f, 0.006512091f,
	0.0069954107f, 0.007499032f, 0.008023193f, 0.008568126f,
	0.009134059f, 0.009721218f, 0.010329823f, 0.010960095f,
	0.011612245f, 0.012286489f, 0.0129830325f, 0.013702083f,
	0.014443845f, 0.015208516f, 0.015996294f, 0.016807377f,
	0.017641956f, 0.018500222f, 0.019382363f, 0.020288564f,
	0.021219011f, 0.022173885f, 0.023153368f, 0.024157634f,
	0.025186861f, 0.026241222f, 0.027320893f, 0.02842604f,
	0.029556835f, 0.030713445f, 0.031896032f, 0.033104766f,
	0.034339808f, 0.035601314f, 0.036889452f, 0.038204372f,
	0.039546236f, 0.0409152f, 0.04231141f, 0.04373503f,
	0.045186203f, 0.046665087f, 0.048171826f, 0.049706567f,
	0.051269464f, 0.05286065f, 0.05448028f, 0.056128494f,
	0.057805438f, 0.059511244f, 0.06124606f, 0.06301002f,
	0.06480327f, 0.066625945f, 0.068478175f, 0.0703601f,
	0.07227185f, 0.07421357f, 0.07618539f, 0.07818743f,
	0.08021983f, 0.082282715f, 0.084376216f, 0.086500466f,
	0.08865559f, 0.09084172f, 0.093058966f, 0.09530747f,
	0.097587354f, 0.09989873f, 0.10224174f, 0.10461649f,
	0.107023105f, 0.10946172f, 0.111932434f, 0.11443538f,
	0.11697067f, 0.119538434f, 0.122138776f, 0.12477182f,
	0.12743768f, 0.13013647f, 0.13286832f, 0.13563333f,
	0.13843162f, 0.14126329f, 0.14412847f, 0.14702727f,
	0.14995979f, 0.15292616f, 0.15592647f, 0.15896083f,
	0.16202939f, 0.1651322f, 0.1682694f, 0.17144111f,
	0.1746474f, 0.17788842f, 0.18116425f, 0.18447499f,
	0.18782078f, 0.19120169f, 0.19461784f, 0.19806932f,
	0.20155625f, 0.20507874f, 0.20863687f, 0.21223076f,
	0.21586053f, 0.21952623f, 0.22322798f, 0.2269659f,
	0.23074007f, 0.23455061f, 0.2383976f, 0.24228115f,
	0.24620135f, 0.2501583f, 0.25415212f, 0.25818288f,
	0.2622507f, 0.26635563f, 0.27049783f, 0.27467734f,
	0.2788943f, 0.28314877f, 0.28744087f, 0.29177067f,
	0.2961383f, 0.3005438f, 0.30498734f, 0.30946895f,
	0.31398875f, 0.3185468f, 0.32314324f, 0.32777813f,
	0.33245155f, 0.33716366f, 0.34191445f, 0.3467041f,
	0.35153264f, 0.35640016f, 0.36130682f, 0.36625263f,
	0.3712377f, 0.37626216f, 0.38132605f, 0.38642946f,
	0.3915725f, 0.39675525f, 0.4019778f, 0.40724024f,
	0.41254264f, 0.4178851f, 0.4232677f, 0.42869052f,
	0.43415368f, 0.4396572f, 0.44520122f, 0.45078582f,
	0.45641103f, 0.46207702f, 0.4677838f, 0.4735315f,
	0.4793202f, 0.48514995f, 0.4910209f, 0.496933f,
	0.5028865f, 0.50888133f, 0.5149177f, 0.5209956f,
	0.52711517f, 0.53327644f, 0.5394795f, 0.5457245f,
	0.55201143f, 0.55834043f, 0.5647115f, 0.57112485f,
	0.57758045f, 0.58407843f, 0.59061885f, 0.5972018f,
	0.60382736f, 0.61049557f, 0.6172066f, 0.62396044f,
	0.63075715f, 0.6375969f, 0.6444797f, 0.65140563f,
	0.65837485f, 0.66538733f, 0.67244315f, 0.6795425f,
	0.6866853f, 0.6938718f, 0.7011019f, 0.7083758f,
	0.71569353f, 0.7230551f, 0.73046076f, 0.73791045f,
	0.74540424f, 0.7529422f, 0.7605245f, 0.76815116f,
	0.7758222f, 0.7835378f, 0.791298f, 0.7991027f,
	0.8069523f, 0.8148466f, 0.82278574f, 0.8307699f,
	0.838799f, 0.8468732f, 0.8549926f, 0.8631572f,
	0.8713671f, 0.8796224f, 0.8879231f, 0.8962694f,
	0.9046612f, 0.91309863f, 0.92158186f, 0.9301109f,
	0.9386857f, 0.9473065f, 0.9559733f, 0.9646863f,
	0.9734453f, 0.9822506f, 0.9911021f, 1.0f,
};

inline int linearToSrgb8bit(float x) {
	if (x <= 0.0f)
		return 0;
	if (x >= 1.0f)
		return 255;
	const float *TABLE = SRGB_8BIT_TO_LINEAR_FLOAT;
	int y = 0;
	for (int i = 128; i != 0; i >>= 1) {
		if (TABLE[y + i] <= x)
			y += i;
	}
	if (x - TABLE[y] <= TABLE[y + 1] - x)
		return y;
	else
		return y + 1;
}


} // namespace

bool resize_bilinear(const float *src, const int32_t src_width,
                     const int32_t src_width_stride, const int32_t src_height,
                     const int32_t channels, const int32_t dst_width,
                     const int32_t dst_width_stride, const int32_t dst_height,
                     std::vector<float> *dst) {
  if (!src) {
    return false;
  }

  if ((src_width < 1) || (src_width_stride < 1) || (src_height < 1)) {
    return false;
  }

  // Up to RGBA
  if ((channels < 1) || (channels > 4)) {
    return false;
  }

  if ((dst_width < 1) || (dst_width_stride < 1) || (dst_height < 1)) {
    return false;
  }

  dst->resize(size_t(dst_width_stride * dst_height * channels));

  int n = stbir_resize_float(
      src, int(src_width), int(src_height), int(src_width_stride), dst->data(),
      int(dst_width), int(dst_height), int(dst_width_stride), int(channels));

  if (n < 1) {
    return false;
  }

  return true;
}

bool imread(const std::string &filename, std::vector<float> *image, int32_t *width, int32_t *height, int32_t *channels, const bool srgb_to_linear)
{
  if (!image) {
    return false;
  }

  int image_width;
  int image_height;
  int n;

   unsigned char *data = stbi_load(filename.c_str(), &image_width, &image_height,
                                    &n, STBI_default);

  if (!data) {
    return false;
  }

  image->resize(size_t(image_width * image_height * n));

  const float inv_scale = 1.0f / 255.0f;

  for (size_t i = 0; i < size_t(image_width * image_height); i++) {
    for (size_t c = 0; c < size_t(n); c++) {

      float value = data[size_t(n) * i + c] * inv_scale;

      // do not apply degamma for alpha channel.
      if (srgb_to_linear && (c < 4)) {
        value = sRGBToLinearFast(value);
      }
      (*image)[size_t(n) * i + c] = value;
    }
  }

  stbi_image_free(data);

  if (width) {
    (*width) = image_width;
  }

  if (height) {
    (*height) = image_height;
  }

  if (channels) {
    (*channels) = n;
  }

  return true;
}

bool imsave(const std::string &filename, std::vector<float> &image, const int32_t width, const int32_t height, const int32_t channels, const bool linear_to_srgb) {

  std::string ext = GetFileExtension(filename);

  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) -> unsigned char {
    return static_cast<unsigned char>(std::tolower(c));
  });

  std::vector<uint8_t> buf(image.size());

  for (size_t i = 0; i < size_t(width * height); i++) {
    for (size_t c = 0; c < size_t(channels); c++) {

      const float src = image[size_t(channels) * i + c];

      int value;

      // do not apply sRGB conversion for alpha channel.
      if (linear_to_srgb && (c < 4)) {
        value = linearToSrgb8bit(src);
      } else {
        value = int(src * 255.0f);
      }

      buf[size_t(channels) * i + c] = static_cast<uint8_t>(std::max(0, std::min(255, value)));
    }
  }

  if ((ext.compare("jpg") == 0) || (ext.compare("jpeg") == 0)) {
    int ret = stbi_write_jpg(filename.c_str(), width, height, channels, reinterpret_cast<void *>(buf.data()), /* quality */100);

    if (ret < 1) {
      return false;
    }
  } else { // save as PNG
    int ret = stbi_write_png(filename.c_str(), width, height, channels, reinterpret_cast<void *>(buf.data()), /* stride */0);

    if (ret < 1) {
      return false;
    }
  }

  return true;

}


}  // namespace nanosnap
