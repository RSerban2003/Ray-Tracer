#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <bvh.h>
DISABLE_WARNINGS_POP()

void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    position = light.endpoint0 + sample * (light.endpoint1 - light.endpoint0);
    color = light.color0 + sample * (light.color1 - light.color0);
}

void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    position = light.v0 + sample.x * light.edge01 + sample.y * light.edge02;
    color = light.color0 + sample.x * (light.color1 - light.color0) + sample.y *
    ((light.color3 + sample.x * (light.color2 - light.color3)) -
    (light.color0 + sample.x * (light.color1 - light.color0)));
}

bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    }
    else {
        //create a ray from the intersection to the light source
        glm::vec3 intersectionPoint = ray.origin + ray.direction * ray.t;
        glm::vec3 shadowDir = glm::normalize(lightPosition - intersectionPoint);
        Ray shadowRay;

        //small offset for the ray so that it doesn't intersect the object it starts from
        shadowRay.origin = intersectionPoint + 0.001f * shadowDir;
        shadowRay.direction = shadowDir;
        shadowRay.t = glm::length(lightPosition - intersectionPoint - 0.001f);

        //check whether the ray intersects anything
        HitInfo collisionHitInfo;
        if (intersectRayWithBVH(state, state.bvh, shadowRay, collisionHitInfo)) {
            //light is not visible
            return false;
        }
        else {
            //light is visible
            return true;
        }
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: implement this function; currently, the light simply passes through
    return lightColor;
}

glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;

    glm::vec3 finalColor = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);
    return computeShading(state, v, l, finalColor, hitInfo);
}

glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    glm::vec3 totalLight(0.0f);

    for (uint32_t i = 0; i < numSamples; ++i)
    {
        //sample the segment light
        glm::vec3 samplePosition, sampleColor;
        sampleSegmentLight(state.sampler.next_1d(), light, samplePosition, sampleColor);

        //create a pointLight object for the pointLight contribution method
        PointLight pointLight;
        pointLight.position = samplePosition;
        pointLight.color = sampleColor;

        //add the contribution of each sample
        totalLight += computeContributionPointLight(state, pointLight, ray, hitInfo);
    }
    return glm::vec3(totalLight.x / numSamples, totalLight.y / numSamples, totalLight.z / numSamples);
}

glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    glm::vec3 totalLight(0.0f);

    for (uint32_t i = 0; i < numSamples; ++i)
    {
        //sample the parallelogram light
        glm::vec3 samplePosition, sampleColor;
        sampleParallelogramLight(state.sampler.next_2d(), light, samplePosition, sampleColor);

        //create a pointLight object for the pointLight contribution method
        PointLight pointLight;
        pointLight.position = samplePosition;
        pointLight.color = sampleColor;

        //add the contribution of each sample
        totalLight += computeContributionPointLight(state, pointLight, ray, hitInfo);
    }
    return glm::vec3(totalLight.x / numSamples, totalLight.y / numSamples, totalLight.z / numSamples);
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}