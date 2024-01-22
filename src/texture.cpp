#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Advanced Task
  // Implement mipmap for trilinear filtering

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation
    u = u > 1.0f ? 1.0f : u;
    u = u < 0.0f ? 0.0f : u;
    v = v > 1.0f ? 1.0f : v;
    v = v < 0.0f ? 0.0f : v;
    auto mipmap = tex.mipmap[level];
    Color c;
    float temp[4];
    int U = round(u * mipmap.width);
    U = U > mipmap.width - 1 ? mipmap.width - 1 : U;
    int V = round(v * mipmap.height);
    V = V > mipmap.height - 1 ? mipmap.height - 1 : V;

    uint8_to_float(temp, &mipmap.texels[4 * (U + V * mipmap.width)]);
    c.r = temp[0];
    c.g = temp[1];
    c.b = temp[2];
    c.a = temp[3];
    return c;
    // return magenta for invalid level
    //return Color(1,0,1,1);

}

static float* lerp(float* a, float* b, float t)
{
    float* c = new float[4];
	c[0] = a[0] * (1.0f - t) + b[0] * t;
	c[1] = a[1] * (1.0f - t) + b[1] * t;
	c[2] = a[2] * (1.0f - t) + b[2] * t;
	c[3] = a[3] * (1.0f - t) + b[3] * t;
	return c;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 4: Implement bilinear filtering
  // return magenta for invalid level
  if (level != 0 || tex.mipmap.size() == 0) return Color(1,0,1,1);

    u = u > 1.0f ? 1.0f : u;
    u = u < 0.0f ? 0.0f : u;
    v = v > 1.0f ? 1.0f : v;
    v = v < 0.0f ? 0.0f : v;
    auto mipmap = tex.mipmap[level];
    
    float U = (u * mipmap.width);
    float V = (v * mipmap.height);
    U = U > mipmap.width - 1 ? mipmap.width - 1 : U;
    V = V > mipmap.height - 1 ? mipmap.height - 1 : V;

    float U_frac = U - floor(U);
    float V_frac = V - floor(V);
    float U_pair = U_frac < 0.5 ? floor(U) - 1 : floor(U) + 1;
    float V_pair = V_frac < 0.5 ? floor(V) - 1 : floor(V) + 1;
    U_pair = U_pair > mipmap.width - 1 ? mipmap.width - 1 : U_pair;
    V_pair = V_pair > mipmap.height - 1 ? mipmap.height - 1 : V_pair;
	U_pair = U_pair < 0 ? 0 : U_pair;
    V_pair = V_pair < 0 ? 0 : V_pair;

    float *temp = new float[4];
    float *temp2 = new float[4];
    float *temp3 = new float[4];
    float *temp4 = new float[4];

    uint8_to_float(temp, &mipmap.texels[4 * (floor(U) + floor(V) * mipmap.width)]);
    uint8_to_float(temp2, &mipmap.texels[4 * (U_pair + floor(V) * mipmap.width)]);
    uint8_to_float(temp3, &mipmap.texels[4 * (floor(U) + V_pair * mipmap.width)]);
    uint8_to_float(temp4, &mipmap.texels[4 * (U_pair + V_pair * mipmap.width)]);

    float t1, t2;
    float sampleU = floor(U) + 0.5f;
    float sampleV = floor(V) + 0.5f;
    t1 = (U - sampleU) / ((U_pair+0.5f) - sampleU);
    t2 = (V - sampleV) / ((V_pair+0.5f) - sampleV);

    auto col1 = lerp(temp2, temp, t1);
    auto col2 = lerp(temp4, temp3, t1);
    auto colfinal = lerp(col1, col1, t2);
    Color c;
    c.r = colfinal[0];
    c.g = colfinal[1];
    c.b = colfinal[2];
    c.a = colfinal[3];
    return c;

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Advanced Task
  // Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248
