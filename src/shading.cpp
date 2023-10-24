#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <shading.h>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded linear gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            { 0.1f, glm::vec3(215.f / 256.f, 210.f / 256.f, 203.f / 256.f) },
            { 0.22f, glm::vec3(250.f / 256.f, 250.f / 256.f, 240.f / 256.f) },
            { 0.5f, glm::vec3(145.f / 256.f, 170.f / 256.f, 175.f / 256.f) },
            { 0.78f, glm::vec3(255.f / 256.f, 250.f / 256.f, 205.f / 256.f) },
            { 0.9f, glm::vec3(170.f / 256.f, 170.f / 256.f, 170.f / 256.f) },
        }
    };

    if (state.features.enableShading) {
        switch (state.features.shadingModel) {
        case ShadingModel::Lambertian:
            return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
        case ShadingModel::Phong:
            return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
        case ShadingModel::BlinnPhong:
            return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
        case ShadingModel::LinearGradient:
            return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Implement basic diffuse shading if you wish to use it
    return sampleMaterialKd(state, hitInfo);
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // TODO: Implement phong shading
    glm::vec3 lightVector = glm::normalize(lightDirection);
    glm::vec3 cameraVector = glm::normalize(cameraDirection);
    glm::vec3 normal = glm::normalize(hitInfo.normal);
    float cosDiffuse = glm::max(glm::dot(lightVector, normal), 0.0f);
    glm::vec3 reflectionVector = 2.0f * glm::dot(lightVector, normal) * normal - lightVector;
    float cosSpecular = glm::max(glm::dot(glm::normalize(reflectionVector), cameraVector), 0.0f);

    glm::vec3 diffuseTerm = lightColor * sampleMaterialKd(state, hitInfo) * cosDiffuse;
    glm::vec3 specularTerm = lightColor * hitInfo.material.ks * glm::pow(cosSpecular, hitInfo.material.shininess);

    return diffuseTerm + specularTerm;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // TODO: Implement blinn-phong shading
    glm::vec3 lightVector = glm::normalize(lightDirection);
    glm::vec3 cameraVector = glm::normalize(cameraDirection);
    glm::vec3 normal = glm::normalize(hitInfo.normal);
    glm::vec3 halfVector = glm::normalize(lightVector + cameraVector);
    float cosDiffuse = glm::max(glm::dot(lightVector, normal), 0.0f);
    float cosSpecular = glm::max(glm::dot(halfVector, normal), 0.0f);

    glm::vec3 diffuseTerm = lightColor * sampleMaterialKd(state, hitInfo) * cosDiffuse;
    glm::vec3 specularTerm = lightColor * hitInfo.material.ks * glm::pow(cosSpecular, hitInfo.material.shininess);

    return diffuseTerm + specularTerm;
}

std::vector<LinearGradient::Component> sort(std::vector<LinearGradient::Component> components)
{
    for (int i = 0; i < components.size(); ++i) {
        float smallestValue = components[i].t;
        int smallestIndex = i;
        int currentIndex = i;
        while (currentIndex < components.size()) {
            float currentValue = components[currentIndex].t;
            if (currentValue < smallestValue) {
                smallestValue = currentValue;
                smallestIndex = currentIndex;
            }
            currentIndex++;
        }
        std::swap(components[i], components[smallestIndex]);
    }
    return components;
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    std::vector<Component> componentsSorted = sort(components);
    if (ti < componentsSorted.front().t)
        return componentsSorted.front().color;
    if (ti > componentsSorted.back().t)
        return componentsSorted.back().color;
    for (Component c : componentsSorted) {
        if (ti == c.t)
            return c.color;
    }

    auto it = std::adjacent_find(components.begin(), components.end(), [ti](const Component& a, const Component& b) {
        return a.t <= ti && ti <= b.t;
    }); // pointer to pair of components

    float leftBoundary = it->t;
    glm::vec3 colorLeft = it->color;
    float rightBoundary = std::next(it)->t;
    glm::vec3 colorRight = std::next(it)->color;

    float weightRight = (ti - leftBoundary) / (rightBoundary - leftBoundary);
    float weightLeft = 1 - weightRight;

    return colorLeft * weightLeft + colorRight * weightRight;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    glm::vec3 lightVector = glm::normalize(lightDirection);
    glm::vec3 normal = glm::normalize(hitInfo.normal);
    float cos_theta = glm::dot(lightVector, normal);
    glm::vec3 kD = gradient.sample(cos_theta);

    cos_theta = glm::max(cos_theta, 0.0f);

    return lightColor * kD * cos_theta;
}