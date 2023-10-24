#include "bvh.h"
#include "draw.h"
#include "extra.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>
#include <stack>

// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex);

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    float maxX = glm::max(glm::max(primitive.v0.position.x, primitive.v1.position.x), primitive.v2.position.x);
    float minX = glm::min(glm::min(primitive.v0.position.x, primitive.v1.position.x), primitive.v2.position.x);
    float maxY = glm::max(glm::max(primitive.v0.position.y, primitive.v1.position.y), primitive.v2.position.y);
    float minY = glm::min(glm::min(primitive.v0.position.y, primitive.v1.position.y), primitive.v2.position.y);
    float maxZ = glm::max(glm::max(primitive.v0.position.z, primitive.v1.position.z), primitive.v2.position.z);
    float minZ = glm::min(glm::min(primitive.v0.position.z, primitive.v1.position.z), primitive.v2.position.z);
    return { .lower = glm::vec3(minX, minY, minZ), .upper = glm::vec3(maxX, maxY, maxZ) };
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    if (primitives.empty())
        return { .lower = glm::vec3(0), .upper = glm::vec3(0) };

    AxisAlignedBox current = computePrimitiveAABB(primitives[0]);
    glm::vec3 lower = current.lower;
    glm::vec3 upper = current.upper;

    for (BVHInterface::Primitive currentPrimitive : primitives) {
        current = computePrimitiveAABB(currentPrimitive);
        lower.x = glm::min(lower.x, current.lower.x);
        upper.x = glm::max(upper.x, current.upper.x);
        lower.y = glm::min(lower.y, current.lower.y);
        upper.y = glm::max(upper.y, current.upper.y);
        lower.z = glm::min(lower.z, current.lower.z);
        upper.z = glm::max(upper.z, current.upper.z);
    }

    return { .lower = lower, .upper = upper };
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    return (primitive.v0.position + primitive.v1.position + primitive.v2.position) / 3.0f;
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    float lenX = aabb.upper.x - aabb.lower.x;
    float lenY = aabb.upper.y - aabb.lower.y;
    float lenZ = aabb.upper.z - aabb.lower.z;

    if (lenX >= lenY && lenX >= lenZ)
        return 0;
    else if (lenY >= lenX && lenY >= lenZ)
        return 1;
    else
        return 2;
}

// forward declaration for helper methods
std::span<BVHInterface::Primitive> sortByAxis(uint32_t axis, std::span<BVHInterface::Primitive> primitives);
float returnAxisValue(uint32_t axis, glm::vec3 vertex);

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    using Primitive = BVHInterface::Primitive;

    primitives = sortByAxis(axis, primitives);
    if (primitives.size() % 2 == 0) {
        return primitives.size() / 2;
    } else {
        return primitives.size() / 2 + 1;
    }
}

// helper (sorting algorithm)
std::span<BVHInterface::Primitive> sortByAxis(uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    for (int i = 0; i < primitives.size(); ++i) {
        float smallestValue = returnAxisValue(axis, computePrimitiveCentroid(primitives[i]));
        int smallestIndex = i;
        int currentIndex = i;
        while (currentIndex < primitives.size()) {
            float currentValue = returnAxisValue(axis, computePrimitiveCentroid(primitives[currentIndex]));
            if (currentValue < smallestValue) {
                smallestValue = currentValue;
                smallestIndex = currentIndex;
            }
            currentIndex++;
        }
        std::swap(primitives[i], primitives[smallestIndex]);
    }

    return primitives;
}

// helper (value for the given axis)
float returnAxisValue(uint32_t axis, glm::vec3 vertex)
{
    if (axis == 0)
        return vertex.x;
    else if (axis == 1)
        return vertex.y;
    else
        return vertex.z;
}

bool intersectRayWithAABB(const AxisAlignedBox& aabb, const Ray& ray)
{
    float tMin = -std::numeric_limits<float>::infinity();
    float tMax = std::numeric_limits<float>::infinity();

    for (int i = 0; i < 3; i++) {
        if (returnAxisValue(i, ray.direction) == 0.0f) {
            // ray is parallel
            if (returnAxisValue(i, ray.origin) < returnAxisValue(i, aabb.lower)
                || returnAxisValue(i, ray.origin) > returnAxisValue(i, aabb.upper))
                return false;
        } else {
            float closerIntersection = (returnAxisValue(i, aabb.lower)
                                           - returnAxisValue(i, ray.origin))
                / returnAxisValue(i, ray.direction);
            float furtherIntersection = (returnAxisValue(i, aabb.upper)
                                            - returnAxisValue(i, ray.origin))
                / returnAxisValue(i, ray.direction);

            if (closerIntersection > furtherIntersection)
                std::swap(closerIntersection, furtherIntersection);

            tMin = std::max(tMin, closerIntersection);
            tMax = std::min(tMax, furtherIntersection);

            if (tMin > tMax)
                return false;
        }
    }
    return true;
}

// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
        // TODO: implement here your (probably stack-based) BVH traversal.
        //
        // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
        // data is not easily extracted. Helper methods are available, however:
        // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
        // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
        // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
        //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
        //
        // In short, you will have to step down the bvh, node by node, and intersect your ray
        // with the node's AABB. If this intersection passes, you should:
        // - if the node is a leaf, intersect with the leaf's primitives
        // - if the node is not a leaf, test the left and right children as well!
        //
        // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
        // and it is likewise possible for a ray to hit both children of a node.
        std::stack<uint32_t> stack;
        stack.push(0u);
        while (!stack.empty()) {
            BVHInterface::Node current = nodes[stack.top()];
            stack.pop();
            if (intersectRayWithAABB(current.aabb, ray)) {
                if (current.isLeaf()) {
                    for (int i = current.primitiveOffset(); i < current.primitiveCount() + current.primitiveOffset(); ++i) {
                        const auto& prim = primitives[i];
                        const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
                        if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                            updateHitInfo(state, prim, ray, hitInfo);
                            is_hit = true;
                        }
                    }
                } else {
                    stack.push(current.leftChild());
                    stack.push(current.rightChild());
                }
            }
        }
    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);

    return is_hit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    Node node;
    // TODO fill in the leaf's data; refer to `bvh_interface.h` for details
    node.aabb = aabb;
    node.data[0] = (node.LeafBit | static_cast<uint32_t>(m_primitives.size()));
    node.data[1] = static_cast<uint32_t>(primitives.size());

    // Copy the current set of primitives to the back of the primitives vector
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));

    return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    Node node;
    // TODO fill in the node's data; refer to `bvh_interface.h` for details
    node.aabb = aabb;
    uint32_t nonLeaf = 0u;
    node.data[0] = (nonLeaf | leftChildIndex);
    node.data[1] = rightChildIndex;
    return node;
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    AxisAlignedBox aabb = computeSpanAABB(primitives);

    // As a starting point, we provide an implementation which creates a single leaf, and stores
    // all triangles inside it. You should remove or comment this, and work on your own recursive
    // construction algorithm that implements the following steps. Make sure to reuse the methods
    // you have previously implemented to simplify this process.
    //
    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    //    (hint; use the `LeafSize` constant)
    // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
    // 3. If it is a node:
    //    3a. Split the range of triangles along the longest axis into left and right subspans,
    //        using either median or SAH-Binning based on the `Features` object
    //    3b. Allocate left/right child nodes
    //        (hint: use `nextNodeIdx()`)
    //    3c. Fill in the current node's data; aabb, left/right child indices
    //    3d. Recursively build left/right child nodes over their respective triangles
    //        (hint; use `std::span::subspan()` to split into left/right ranges)

    // Just configure the current node as a giant leaf for now

    if (primitives.size() <= this->LeafSize) {
        // it is a leaf
        this->m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
    } else {
        // it is not a leaf
        uint32_t longestAxis = computeAABBLongestAxis(aabb);
        size_t splitIndex = splitPrimitivesByMedian(aabb, longestAxis, primitives);

        uint32_t leftChildIndex = nextNodeIdx();
        uint32_t rightChildIndex = nextNodeIdx();

        this->m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftChildIndex, rightChildIndex);

        buildRecursive(scene, features, primitives.subspan(0, splitIndex), leftChildIndex);
        buildRecursive(scene, features, primitives.subspan(splitIndex), rightChildIndex);
    }
}

uint32_t maxDepth(BVHInterface::Node node, std::vector<BVHInterface::Node> nodes);

// TODO: Standard feature, or part of it
// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
    m_numLevels = maxDepth(m_nodes[RootIndex], m_nodes);
}

uint32_t maxDepth(BVHInterface::Node node, std::vector<BVHInterface::Node> nodes)
{
    if (node.isLeaf())
        return 1;

    uint32_t depthL = maxDepth(nodes[node.leftChild()], nodes);
    uint32_t depthR = maxDepth(nodes[node.rightChild()], nodes);

    if (depthL > depthR)
        return depthL + 1;
    else
        return depthR + 1;
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    uint32_t counter = 0u;
    std::stack<uint32_t> stack;
    stack.push(0);
    while (!stack.empty()) {
        BVHInterface::Node current = m_nodes[stack.top()];
        stack.pop();

        if (current.isLeaf()) {
            counter += static_cast<uint32_t>(1u);
            continue;
        }

        stack.push(current.leftChild());
        stack.push(current.rightChild());
    }

    m_numLeaves = counter;
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
    if (level >= m_numLevels) {
        return;
    }

    std::stack<std::tuple<uint32_t, int>> nodeStack;
    nodeStack.push(std::tuple(RootIndex, 0));

    while (!nodeStack.empty()) {
        std::tuple currTuple = nodeStack.top();
        nodeStack.pop();
        Node current = m_nodes[std::get<0>(currTuple)];
        if (std::get<1>(currTuple) == level) {
            drawAABB(current.aabb, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.5f);
        }
        if (!current.isLeaf()) {
            nodeStack.push(std::tuple(current.leftChild(), std::get<1>(currTuple) + 1));
            nodeStack.push(std::tuple(current.rightChild(), std::get<1>(currTuple) + 1));
        }
    }
}

void BVH::debugDrawLeaf(int leafIndex)
{
    if (leafIndex >= m_numLeaves) {
        return;
    }

    uint32_t startIndexLeaf = m_nodes[RootIndex].primitiveOffset();
    for (uint32_t i = 0; i < leafIndex; ++i) {
        startIndexLeaf += m_nodes[RootIndex].primitiveCount();
    }
    uint32_t endIndexLeaf = startIndexLeaf + m_nodes[RootIndex].primitiveCount();

    for (uint32_t i = startIndexLeaf; i < endIndexLeaf; ++i) {
        Primitive primitive = m_primitives[i];
        drawTriangle(primitive.v0, primitive.v1, primitive.v2);
    }
}
