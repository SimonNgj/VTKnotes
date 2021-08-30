// The implementation of MeshLab Filters/Sampling/Poisson-disk sampling function
// Load a mesh and then sample it to get a point cloud


#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export_ply.h>

#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
    Use<MyEdge>     ::AsEdgeType,
    Use<MyFace>     ::AsFaceType> {};

class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  > {};
class MyFace : public Face< MyUsedTypes, face::FFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge : public Edge<MyUsedTypes> {};
class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace>, vector<MyEdge>  > {};

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        printf("Usage trimesh_base <meshfilename.obj> radius\n");
        return -1;
    }

    MyMesh m;

    if (tri::io::Importer<MyMesh>::Open(m, argv[1]) != 0)
    {
        printf("Error reading file  %s\n", argv[1]);
        exit(0);
    }
    tri::SurfaceSampling<MyMesh, tri::TrivialSampler<MyMesh> >::SamplingRandomGenerator().initialize(time(0));

    //----------------------------------------------------------------------
    // Advanced Sample
    // Make a feature dependent Poisson Disk sampling
    MyMesh MontecarloSurfaceMesh;

    std::vector<Point3f> sampleVec;
    tri::TrivialSampler<MyMesh> mps(sampleVec);
    tri::SurfaceSampling<MyMesh, tri::TrivialSampler<MyMesh> >::Montecarlo(m, mps, 10000);
    tri::BuildMeshFromCoordVector(MontecarloSurfaceMesh, sampleVec);
    tri::io::ExporterPLY<MyMesh>::Save(MontecarloSurfaceMesh, "RandomMesh.ply");

    return 0;
}
