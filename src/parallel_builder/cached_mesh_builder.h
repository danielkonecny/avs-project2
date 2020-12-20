/**
 * @file    cached_mesh_builder.h
 *
 * @author  Daniel Konecny <xkonec75@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using pre-computed field
 *
 * @date    20 December 2020
 **/

#ifndef CACHED_MESH_BUILDER_H
#define CACHED_MESH_BUILDER_H

#include <vector>
#include <limits>
#include "base_mesh_builder.h"

/**
 * @brief The CachedMeshBuilder class
 */
class CachedMeshBuilder : public BaseMeshBuilder
{
public:
    CachedMeshBuilder(unsigned gridEdgeSize);

protected:
	float computeDistance(const Vec3_t<float> &pos, const ParametricScalarField &field);
	void computeCachedDistances(const ParametricScalarField &field);
    unsigned marchCubes(const ParametricScalarField &field);
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data(); }

    std::vector<Triangle_t> mTriangles; ///< Temporary array of triangles
    std::vector<float> cachedCubes;
};

#endif // CACHED_MESH_BUILDER_H
