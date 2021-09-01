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

	std::string namefile = "../files/testVCGclean";
	string meshFilename = namefile + ".ply";

	bool openFile2 = tri::io::Importer<MyMesh>::Open(original, meshFilename.c_str());

	////////////////// Remove Isolated pieces (wrt Face Num.) /////////////////////
	vcg::face::FFAdj<MyMesh>::FFAdj();
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(original);

	int minCC_f = 50;
	std::pair<int, int> delInfo_f = tri::Clean<MyMesh>::RemoveSmallConnectedComponentsSize(original, minCC_f);

	//"removeUnref"
	int delvert_noFace = tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
	tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
	Allocator<MyMesh>::CompactEveryVector(original);
	std::string meshFilename2a = namefile + "_removedFaces_noFaces.ply";
	vcg::tri::io::Exporter<MyMesh>::Save(original, meshFilename2a.c_str());

 //   /////////////////////// Remove Isolated pieces (diameter) ////////////////////////
 //   float minCC_d = 5;
 //   std::pair<int, int> delInfo_d = tri::Clean<MyMesh>::RemoveSmallConnectedComponentsDiameter(original, minCC_d);
 //   int delvert_diameter = tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
 //   tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
 //   Allocator<MyMesh>::CompactEveryVector(original);
 //
	//std::string meshFilename2b = namefile + "_removedFaces_diameter.ply";
	//vcg::tri::io::Exporter<MyMesh>::Save(original, meshFilename2b.c_str());
	//////////////////////////////////////////

	return 0;
}