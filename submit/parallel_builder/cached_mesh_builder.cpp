/**
 * @file    cached_mesh_builder.cpp
 *
 * @author  Daniel Konecny <xkonec75@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using pre-computed field
 *
 * @date    18 December 2020
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "cached_mesh_builder.h"

CachedMeshBuilder::CachedMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Cached")
{

}

float CachedMeshBuilder::computeDistance(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
	// NOTE: This method is called from "buildCube(...)"!

    // 1. Store pointer to and number of 3D points in the field
    //    (to avoid "data()" and "size()" call in the loop).
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    for(unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void CachedMeshBuilder::computeCachedDistances(const ParametricScalarField &field)
{
	size_t cacheSize = (mGridSize+1)*(mGridSize+1)*(mGridSize+1);
	cachedCubes = std::vector<float>(cacheSize, 0.0);

	#pragma omp parallel for schedule(guided)
	for(size_t i = 0; i < cacheSize; i++)
	{
		float distance;
		size_t fieldX = (size_t) (i % (mGridSize+1));
	    size_t fieldY = (size_t) ((i / (mGridSize+1)) % (mGridSize+1));
	    size_t fieldZ = (size_t) (i / ((mGridSize+1) * (mGridSize+1)));
	    
		Vec3_t<float> position(fieldX*mGridResolution,
							   fieldY*mGridResolution,
							   fieldZ*mGridResolution);

		distance = computeDistance(position, field);
		
		#pragma omp critical
		cachedCubes[i] = distance;
	}
}

unsigned CachedMeshBuilder::marchCubes(const ParametricScalarField &field)
{
	computeCachedDistances(field);

    // 1. Compute total number of cubes in the grid.
    size_t totalCubesCount = mGridSize*mGridSize*mGridSize;

    unsigned totalTriangles = 0;

    // 2. Loop over each coordinate in the 3D grid.
    #pragma omp parallel for reduction(+:totalTriangles) schedule(guided)
    for(size_t i = 0; i < totalCubesCount; ++i)
    {
        // 3. Compute 3D position in the grid.
        Vec3_t<float> cubeOffset( i % mGridSize,
                                 (i / mGridSize) % mGridSize,
                                  i / (mGridSize*mGridSize));

        // 4. Evaluate "Marching Cube" at given position in the grid and
        //    store the number of triangles generated.
        totalTriangles += buildCube(cubeOffset, field);
    }

    // 5. Return total number of triangles generated.
    return totalTriangles;
}

float CachedMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // NOTE: This method is called from "buildCube(...)"!
    size_t fieldX = (size_t) floor(pos.x / mGridResolution + 0.5);
    size_t fieldY = (size_t) floor(pos.y / mGridResolution + 0.5);
    size_t fieldZ = (size_t) floor(pos.z / mGridResolution + 0.5);
    size_t fieldIndex = fieldX + fieldY*(mGridSize+1) + fieldZ*(mGridSize+1)*(mGridSize+1);

    return cachedCubes[fieldIndex];
}

void CachedMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
	// NOTE: This method is called from "buildCube(...)"!

    // Store generated triangle into vector (array) of generated triangles.
    // The pointer to data in this array is return by "getTrianglesArray(...)" call
    // after "marchCubes(...)" call ends.
    #pragma omp critical
    mTriangles.push_back(triangle);
}
