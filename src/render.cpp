#include "render.h"
#include "bvh_interface.h"
#include "draw.h"
#include "extra.h"
#include "light.h"
#include "recursive.h"
#include "sampler.h"
#include "screen.h"
#include "shading.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

// This function is provided as-is. You do not have to implement it.
// Given relevant objects (scene, bvh, camera, etc) and an output screen, multithreaded fills
// each of the pixels using one of the below `renderPixel*()` functions, dependent on scene
// configuration. By default, `renderPixelNaive()` is called.
void renderImage(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    // Either directly render the image, or pass through to extra.h methods
    if (features.extra.enableDepthOfField) {
        renderImageWithDepthOfField(scene, bvh, features, camera, screen);
    } else if (features.extra.enableMotionBlur) {
        renderImageWithMotionBlur(scene, bvh, features, camera, screen);
    } else {
#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
        for (int y = 0; y < screen.resolution().y; y++) {
            for (int x = 0; x != screen.resolution().x; x++) {
                // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
                // Note; we seed the sampler for consistenct behavior across frames
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };
                auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
                auto L = renderRays(state, rays);
                screen.setPixel(x, y, L);
            }
        }
    }

    // Pass through to extra.h for post processing
    if (features.extra.enableBloomEffect) {
        postprocessImageWithBloom(scene, features, camera, screen);
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples for this pixel.
// This method forwards to `generatePixelRaysMultisampled` and `generatePixelRaysStratified` when necessary.
std::vector<Ray> generatePixelRays(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    if (state.features.numPixelSamples > 1) {
        if (state.features.enableJitteredSampling) {
            return generatePixelRaysStratified(state, camera, pixel, screenResolution);
        } else {
            return generatePixelRaysMultisampled(state, camera, pixel, screenResolution);
        }
    } else {
        // Generate single camera ray placed at the pixel's center
        // Note: (-1, -1) at the bottom left of the screen,
        //       (+1, +1) at the top right of the screen.
        glm::vec2 position = (glm::vec2(pixel) + 0.5f) / glm::vec2(screenResolution) * 2.f - 1.f;
        return { camera.generateRay(position) };
    }
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// uniformly throughout this pixel.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is unit-tested, so do not change the function signature.
std::vector<Ray> generatePixelRaysMultisampled(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples camera rays uniformly distributed across the pixel. Use
    auto numSamples = state.features.numPixelSamples;
    std::vector<Ray> rays;

    /*
    I assumed the number of samples does not have to be a perfect
    square, so this approach, even if it might not place the samples
    perfectly uniform, the difference should be as unnoticeable as
    possible and ALL the samples should be placed, even for numbers
    with inconvenient factorizations, such as prime numbers (that would result in poor
    row/column-wise distributions).
    PS: If the number of samples is in fact a perfect square then the samples will be placed perfectly uniform.
    */

    //determine the size of the cell each sample occupies
    float cellSize = 1.0f / std::sqrt(numSamples);

    //offset the samples to the center of the cell, and place the cells column-wise until the width of the pixel is reached, then switch to a new row
    float xOffset, yOffset;
    yOffset = 0.5f * cellSize;
    while (yOffset < 1.0f) {
        xOffset = 0.5f * cellSize;
        while(xOffset < 1.0f && rays.size() < numSamples) {

            //calculate the ray's position on the screen: add the pixel's coordinates, normalize, scale and translate to screen space
            glm::vec2 position = (glm::vec2(pixel) + glm::vec2(xOffset, yOffset)) / glm::vec2(screenResolution) * 2.f - 1.f;

            //add the ray to the vector
            rays.push_back(camera.generateRay(position));
            xOffset += cellSize;
        }

        yOffset += cellSize;
    }
    return rays;
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// using jittered sampling throughout this pixel. Given NxN cells across one pixel, each ray sample is randomly
// placed somewhere within a cell.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
std::vector<Ray> generatePixelRaysStratified(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples * numSamples camera rays as jittered samples across the pixel.
    auto numSamples = static_cast<uint32_t>(std::round(std::sqrt(float(state.features.numPixelSamples))));
    std::vector<Ray> rays;

    //calculate the width of each cell based on the square root of the number of samples
    float cellSize = 1.0f / numSamples;

    //loop through all the cells
    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numSamples; j++) {

            //calculate random x and y offset for a jittered effect and scale it to cell size
            float xOffset = (i + state.sampler.next1d()) * cellSize;
            float yOffset = (j + state.sampler.next1d()) * cellSize;

            //calculate the ray's position on the screen: add the pixel's coordinates, normalize, scale and translate to screen space
            glm::vec2 pos = (glm::vec2(pixel) + glm::vec2(xOffset, yOffset)) / glm::vec2(screenResolution) * 2.f - 1.f;

            //add the ray to the vector
            rays.push_back(camera.generateRay(pos));
        }
    }
    return rays;
}