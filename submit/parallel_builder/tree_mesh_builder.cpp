/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  Daniel Konecny <xkonec75@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    18 December 2020
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"


TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned TreeMeshBuilder::treeDive(const ParametricScalarField &field,
								   Vec3_t<float> position, unsigned depth)
{
	const unsigned treeGridSize = 2;
	size_t treeCubesCount = treeGridSize*treeGridSize*treeGridSize;

	unsigned totalTriangles = 0;

	if(depth >= mGridSize)
	{
		unsigned tmpTriangles = buildCube(position, field);

		#pragma omp atomic update
		totalTriangles += tmpTriangles;
	}
	else
	{
		unsigned edgeLength = mGridSize/depth;
		unsigned halfEdge = edgeLength/2;
		float distance = field.getIsoLevel() + (sqrt(3)/2) * (edgeLength*mGridResolution);

		Vec3_t<float> center((position.x + halfEdge)*mGridResolution,
							 (position.y + halfEdge)*mGridResolution,
							 (position.z + halfEdge)*mGridResolution);

		if(evaluateFieldAt(center, field) <= distance)
		{
			for(size_t i = 0; i < treeCubesCount; i++)
			{
				#pragma omp task shared(totalTriangles)
				{
					Vec3_t<float> pos(position.x + ( i % treeGridSize)*halfEdge,
									  position.y + ((i / treeGridSize) % treeGridSize)*halfEdge,
									  position.z + ( i / (treeGridSize * treeGridSize))*halfEdge);

					unsigned tmpTriangles = treeDive(field, pos, treeGridSize*depth);

					#pragma omp atomic update
					totalTriangles += tmpTriangles;
				}
			}
		}
	}
	#pragma omp taskwait
	return totalTriangles;
}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.
    unsigned totalTriangles = 0;
    Vec3_t<float> origin(0, 0, 0);
    #pragma omp parallel
	{
		#pragma omp single
		{
    		totalTriangles = treeDive(field, origin, 1);
		}
	}
    return totalTriangles;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
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

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    // NOTE: This method is called from "buildCube(...)"!

    // Store generated triangle into vector (array) of generated triangles.
    // The pointer to data in this array is return by "getTrianglesArray(...)" call
    // after "marchCubes(...)" call ends.
    #pragma omp critical
    mTriangles.push_back(triangle);
}
