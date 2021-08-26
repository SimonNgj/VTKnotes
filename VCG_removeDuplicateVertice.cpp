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
	MyMesh original;

	bool loaddata = tri::io::Importer<MyMesh>::Open(original, "../files/testInverseFace1.ply");

	////////////////// Remove duplicated vertices /////////////////////
	int delvert = tri::Clean<MyMesh>::RemoveDuplicateVertex(original);
	if (delvert != 0) {
		tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(original);
		tri::UpdateBounding<MyMesh>::Box(original);
	}

	////////////////// Save ///////////////////////
	std::string meshFilename2 = "../files/testInverseFace_result.ply";
	vcg::tri::io::Exporter<MyMesh>::Save(original, meshFilename2.c_str());

	return 0;
}
