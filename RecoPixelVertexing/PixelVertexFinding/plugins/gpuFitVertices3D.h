#ifndef RecoPixelVertexing_PixelVertexFinding_plugins_gpuFitVertices3D_h
#define RecoPixelVertexing_PixelVertexFinding_plugins_gpuFitVertices3D_h

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "HeterogeneousCore/CUDAUtilities/interface/HistoContainer.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cuda_assert.h"

#include "gpuVertexFinder.h"

namespace gpuVertexFinder {

  __device__ __forceinline__ void fitVertices3D(ZVertices* pdata,
                                              WorkSpace* pws,
                                              float chi2Max  // for outlier rejection
  ) {
    constexpr bool verbose = false;  // in principle the compiler should optmize out if false

    auto& __restrict__ data = *pdata;
    auto& __restrict__ ws = *pws;
    auto nt = ws.ntrks;
    float const* __restrict__ zt = ws.zt;
    float const* __restrict__ ezt2 = ws.ezt2;

    float* __restrict__ zv = data.zv;
    float* __restrict__ wv = data.wv;
    
    float* __restrict__ xv = data.xv;
    float* __restrict__ wvx = data.wvx;

    float* __restrict__ yv = data.yv;
    float* __restrict__ wvy = data.wvy;
    
    float* __restrict__ chi2 = data.chi2;
    uint32_t& nvFinal = data.nvFinal;
    uint32_t& nvIntermediate = ws.nvIntermediate;

    int32_t* __restrict__ nn = data.ndof;
    int32_t* __restrict__ iv = ws.iv;

    assert(pdata);
    assert(zt);
    assert(xt);
    assert(yt);


    assert(nvFinal <= nvIntermediate);
    nvFinal = nvIntermediate;
    auto foundClusters = nvFinal;

    // zero
    for (auto i = threadIdx.x; i < foundClusters; i += blockDim.x) {
      zv[i] = 0;
      wv[i] = 0;
      xv[i] = 0;
      wvx[i] = 0;
      uv[i] = 0;
      wvy[i] = 0;
      chi2[i] = 0;
    }

    // only for test
    __shared__ int noise;
    if (verbose && 0 == threadIdx.x)
      noise = 0;

    __syncthreads();

    // compute cluster location
    for (auto i = threadIdx.x; i < nt; i += blockDim.x) {
      if (iv[i] > 9990) {
        if (verbose)
          atomicAdd(&noise, 1);
        continue;
      }
      assert(iv[i] >= 0);
      assert(iv[i] < int(foundClusters));
      auto w = 1.f / ezt2[i];
      atomicAdd_block(&zv[iv[i]], zt[i] * w);
      atomicAdd_block(&wv[iv[i]], w);
    }

    __syncthreads();
    // reuse nn
    for (auto i = threadIdx.x; i < foundClusters; i += blockDim.x) {
      assert(wv[i] > 0.f);
      zv[i] /= wv[i];
      xv[i] /= wvx[i];
      yv[i] /= wvy[i];
      nn[i] = -1;  // ndof
    }
    __syncthreads();

    // compute chi2
    for (auto i = threadIdx.x; i < nt; i += blockDim.x) {
      if (iv[i] > 9990)
        continue;

      auto c2 = zv[iv[i]] - zt[i];  // output z position of the vertex index - input z position of track 
      c2 *= c2 / ezt2[i]; // ez2 input 
      if (c2 > chi2Max) {
        iv[i] = 9999;
        continue;
      }
      if (iv[i] > 9990)
        continue;

      auto c2 = xv[iv[i]] - xt[i];  // output x position of the vertex index - input x position of track 
      c2 *= c2 / ext2[i]; // ex2 input 
      if (c2 > chi2Max) {
        iv[i] = 9999;
        continue;
      }
      if (iv[i] > 9990)
        continue;

      auto c2 = yv[iv[i]] - yt[i];  // output y position of the vertex index - input y position of track 
      c2 *= c2 / eyt2[i]; // ey2 input 
      if (c2 > chi2Max) {
        iv[i] = 9999;
        continue;
      }
      atomicAdd_block(&chi2[iv[i]], c2);
      atomicAdd_block(&nn[iv[i]], 1);
    }
    __syncthreads();
    for (auto i = threadIdx.x; i < foundClusters; i += blockDim.x)
      if (nn[i] > 0)
        wv[i] *= float(nn[i]) / chi2[i];

    if (verbose && 0 == threadIdx.x)
      printf("found %d proto clusters ", foundClusters);
    if (verbose && 0 == threadIdx.x)
      printf("and %d noise\n", noise);
  }

  __global__ void fitVerticesKernel3D(ZVertices* pdata,
                                    WorkSpace* pws,
                                    float chi2Max  // for outlier rejection
  ) {
    printf("[gpuFitVertices3D] running 3D vertex fitting function"); 
    fitVertices(pdata, pws, chi2Max);
  }

}  // namespace gpuVertexFinder

#endif  // RecoPixelVertexing_PixelVertexFinding_plugins_gpuFitVertices_h
