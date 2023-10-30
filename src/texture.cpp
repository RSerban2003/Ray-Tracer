#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.

// Note: the center of the first pixel should be at coordinates (0.5, 0.5)
// Given texcoords, return the corresponding pixel of the image
// The pixel are stored in a 1D array of row major order
// you can convert from position (i,j) to an index using the method seen in the lecture
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    //convert texture coordinates to image coordinates and add the center offset
    float i = image.width * texCoord.x - 0.5f;
    float j = image.height * texCoord.y - 0.5f;

    //convert the 2D indices to a 1D index using the method from the lecture and round the coordinates to obtain the nearest texel
    int index = std::round(j) * image.width + std::round(i);

    //we return the index and if it is out of bounds we wrap around
    return image.pixels[index % image.pixels.size()];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.

// Given texcoords, return the corresponding pixel of the image
// The pixel are stored in a 1D array of row major order
// you can convert from position (i,j) to an index using the method seen in the lecture
// Note, the center of the first pixel is at image coordinates (0.5, 0.5)
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    //convert texture coordinates to image coordinates
    float i = (texCoord.x * image.width) - 0.5f;
    float j = (texCoord.y * image.height) - 0.5f;

    //determine the integer parts of the coordinates and calculate the interpolation weights
    int leftTexel = std::floor(i);
    int rightTexel = leftTexel + 1;
    int topTexel = std::floor(j);
    int bottomTexel = topTexel + 1;

    float weightLeft = (rightTexel - i);
    float weightRight = 1.0f - weightLeft;
    float weightTop = (bottomTexel - j);
    float weightBottom = 1.0f - weightTop;

    //get the neighbouring texels (in case of out-of-bounds index wrap around the texture)
    glm::vec3 topLeft = image.pixels[(topTexel * image.width + leftTexel) % image.pixels.size()];
    glm::vec3 topRight = image.pixels[(topTexel * image.width + rightTexel) % image.pixels.size()];
    glm::vec3 bottomLeft = image.pixels[(bottomTexel * image.width + leftTexel) % image.pixels.size()];
    glm::vec3 bottomRight = image.pixels[(bottomTexel * image.width + rightTexel) % image.pixels.size()];

    //return the bilinearly interpolated texel (the helper method from the interpolate class only uses 3 coordinates not 4, so it is not helpful here)
    return weightLeft * weightTop * topLeft + weightRight * weightTop * topRight + weightLeft * weightBottom * bottomLeft + weightRight * weightBottom * bottomRight;
}