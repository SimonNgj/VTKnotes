/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file trimesh_normal.cpp
\ingroup code_sample

\brief An example of all the methods for computing normals over a mesh.

*/
#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
	Use<MyEdge>     ::AsEdgeType,
	Use<MyFace>     ::AsFaceType> {};

class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::BitFlags, vertex::Mark> {};
class MyFace : public Face< MyUsedTypes, face::Mark, face::VertexRef, face::VFAdj, face::FFAdj, face::Normal3f, face::BitFlags > {};
class MyEdge : public Edge<MyUsedTypes> {};
class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace>, vector<MyEdge>  > {};

int main(int argc, char** argv)
{
	MyMesh original, toremesh;

	bool loaddata = tri::io::Importer<MyMesh>::Open(original, "../files/testInverseFace1.ply");

	////float targetLenPerc = .2f;
	////int iterNum = 20;
	////float creaseAngle = 30.f;
	////float maxSurfDistPerc = 0.001f;

	//// Mesh cleaning
	//tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
	//Allocator<MyMesh>::CompactEveryVector(original);

	//tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(original);
	//tri::UpdateBounding<MyMesh>::Box(original);

	//vcg::tri::Append<MyMesh, MyMesh>::MeshCopy(toremesh, original);
	//tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(toremesh);
	//tri::UpdateBounding<MyMesh>::Box(toremesh);

	//// Start to flip all normals to outside
	////vcg::face::FFAdj<MyMesh>::FFAdj();
	//tri::UpdateTopology<MyMesh>::FaceFace(toremesh);

	//bool oriented, orientable;
	////if (vcg::tri::Clean<MyMesh>::CountNonManifoldEdgeFF(m) > 0) {
	////    std::cout << "Mesh has some not 2-manifold faces, Orientability requires manifoldness" << std::endl; // text
	////    return 0; // can't continue, mesh can't be processed
	////}
	//vcg::tri::Clean<MyMesh>::OrientCoherentlyMesh(toremesh, oriented, orientable);
	//vcg::tri::Clean<MyMesh>::FlipNormalOutside(toremesh);
	////vcg::tri::Clean<MyMesh>::FlipMesh(toremesh);
	////vcg::tri::UpdateTopology<MyMesh>::FaceFace(toremesh);
	////vcg::tri::UpdateTopology<MyMesh>::TestFaceFace(toremesh);
	//vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(toremesh);
	//vcg::tri::UpdateNormal<MyMesh>::PerVertexFromCurrentFaceNormal(toremesh);

	//////////////////// Remove Isolated pieces (wrt size.) /////////////////////
	//float minCC = 0.2;
	//std::pair<int, int> delInfo = tri::Clean<MyMesh>::RemoveSmallConnectedComponentsDiameter(toremesh, minCC);
	//int delvert = tri::Clean<MyMesh>::RemoveUnreferencedVertex(toremesh);


	////////////////// Remove Isolated pieces (wrt Face Num.) /////////////////////
	//int minCC = 25;
	//std::pair<int, int> delInfo = tri::Clean<MyMesh>::RemoveSmallConnectedComponentsSize(toremesh, minCC);

	////"removeUnref"
	//int delvert = tri::Clean<MyMesh>::RemoveUnreferencedVertex(toremesh);
	//tri::Clean<MyMesh>::RemoveUnreferencedVertex(toremesh);
	//Allocator<MyMesh>::CompactEveryVector(toremesh);

	//////////////////// Remove duplicated face /////////////////////
	//int total_faces = tri::Clean<MyMesh>::RemoveDuplicateFace(toremesh);

	//////////////////// Remove duplicated vertices /////////////////////
	//int delvert = tri::Clean<MyMesh>::RemoveDuplicateVertex(toremesh);
	//if (delvert != 0) {
	//	tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(toremesh);
	//	tri::UpdateBounding<MyMesh>::Box(toremesh);
	//}

	//////////////////// Inverse face /////////////////////
	bool flipped = true;
	bool onlySelected = false;

	if (flipped)
		tri::Clean<MyMesh>::FlipMesh(original, onlySelected);
	else
		flipped = tri::Clean<MyMesh>::FlipNormalOutside(original);

	tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(original);
	tri::UpdateBounding<MyMesh>::Box(original);
	//m.clearDataMask(MeshModel::MM_FACEFACETOPO);

	////////////////// Save ///////////////////////
	std::string meshFilename2 = "../files/testInverseFace_result.ply";
	vcg::tri::io::Exporter<MyMesh>::Save(original, meshFilename2.c_str());



	//float lengthThr = targetLenPerc * (original.bbox.Diag() / 100.f);
	//float maxSurfDist = maxSurfDistPerc * (original.bbox.Diag() / 100.f);
	//printf("Length Thr: %8.3f ~ %4.2f %% on %5.3f\n", lengthThr, targetLenPerc, original.bbox.Diag());

	//IsotropicRemeshing<MyMesh>::Params params;
	//params.SetTargetLen(lengthThr);
	//params.SetFeatureAngleDeg(creaseAngle);
	//params.iter = iterNum;

	//if (maxSurfDistPerc != 0)
	//{
	//	params.surfDistCheck = true;
	//	params.maxSurfDist = maxSurfDist;
	//}
	//else
	//{
	//	params.surfDistCheck = false;
	//}

	//params.cleanFlag = true;
	//params.userSelectedCreases = false;

	//printf(" Input mesh %8i v %8i f\n", toremesh.VN(), toremesh.FN());
	//IsotropicRemeshing<MyMesh>::Do(toremesh, original, params);
	//vcg::tri::io::ExporterPLY<MyMesh>::Save(toremesh, "../files/remesh.ply");
	//printf("Output mesh %8i v %8i f\n", toremesh.VN(), toremesh.FN());

	return 0;
}
