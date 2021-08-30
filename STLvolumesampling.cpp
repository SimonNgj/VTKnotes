/*#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/import_stl.h>
#include<wrap/io_trimesh/export_off.h>
#include<wrap/io_trimesh/export_ply.h>
#include<wrap/io_trimesh/export_dxf.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/voronoi_processing.h>
#include<vcg/complex/algorithms/voronoi_volume_sampling.h>


using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
    Use<MyEdge>     ::AsEdgeType,
    Use<MyFace>     ::AsFaceType> {};

class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  > {};
class MyFace : public Face< MyUsedTypes, face::VertexRef, face::Normal3f, face::BitFlags, face::Mark, face::VFAdj, face::FFAdj > {};
class MyEdge : public Edge< MyUsedTypes, edge::VertexRef, edge::BitFlags> {};
class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyEdge>, vector<MyFace>   > {};

class EmEdge;
class EmFace;
class EmVertex;
struct EmUsedTypes : public UsedTypes<	Use<EmVertex>   ::AsVertexType,
    Use<EmEdge>     ::AsEdgeType,
    Use<EmFace>     ::AsFaceType> {};

class EmVertex : public Vertex<EmUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  > {};
class EmFace : public Face< EmUsedTypes, face::VertexRef, face::BitFlags, face::VFAdj > {};
class EmEdge : public Edge< EmUsedTypes, edge::VertexRef> {};
class EmMesh : public tri::TriMesh< vector<EmVertex>, vector<EmEdge>, vector<EmFace>   > {};


int main(int argc, char** argv)
{
    MyMesh mOff; // the offsetted surface

    tri::io::ImporterPLY<MyMesh>::Open(mOff, "../files/10.STL");
    tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(mOff);

    std::vector<Point3f> sampleSurVec, sampleOffVec;
    MontecarloSampling(mOff, sampleOffVec, 10000);
    sampleSurVec.insert(sampleSurVec.end(), sampleOffVec.begin(), sampleOffVec.end());


    printf("Read %i vn %i fn \n", mOff.vn, mOff.fn);

    float poissonRadius = mOff.bbox.Diag() * atof(argv[3]);
    //  float poissonRadius = mOff.bbox.Diag()/10; 
    printf("Poisson Radius %f\n", poissonRadius);

    float sampleSurfRadius = mOff.bbox.Diag() / 100.0f;
    int montecarloSampleNum = 100000;

    MyMesh  seedM;
    VoronoiVolumeSampling<MyMesh> vvs(mOff, seedM);
    printf("Sampling Surface at a radius %f ", sampleSurfRadius);
    vvs.Init(sampleSurfRadius);
    tri::BuildMeshFromCoordVector(vvs.seedDomainMesh, sampleSurVec);
    printf("Sampled\n");
    vvs.BuildVolumeSampling(montecarloSampleNum, 0, poissonRadius);


    tri::io::ExporterPLY<MyMesh>::Save(vvs.seedDomainMesh, "seedDomainMesh.ply");
    //tri::io::ExporterPLY<MyMesh>::Save(vvs.poissonSurfaceMesh, "poissonSurfaceMesh.ply");
    tri::io::ExporterPLY<MyMesh>::Save(vvs.montecarloVolumeMesh, "montecarloVolumeMesh.ply");
    tri::io::ExporterPLY<MyMesh>::Save(seedM, "seedMesh0.ply");
    //  vvs.restrictedRelaxationFlag=true;
    //  vvs.BarycentricRelaxVoronoiSamples(10);
    vvs.QuadricRelaxVoronoiSamples(10);

    tri::UpdateColor<MyMesh>::PerVertexQualityRamp(seedM);
    tri::io::ExporterPLY<MyMesh>::Save(seedM, "seedMesh1.ply", tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);
    vvs.QuadricRelaxVoronoiSamples(10);

    tri::UpdateColor<MyMesh>::PerVertexQualityRamp(seedM);
    tri::io::ExporterPLY<MyMesh>::Save(seedM, "seedMesh2.ply", tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);
    printf("\n Saved %i points \n", seedM.vn);
    // Second Pipeline 
    return true;

    return 0;
}*/

/*
#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/clustering.h>

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
    if (argc < 3)
    {
        printf("Usage trimesh_base <meshfilename> radius (as perc) (\n");
        return -1;
    }

    MyMesh m;
    MyMesh subM;
    MyMesh cluM;
    MyMesh rndM;

    tri::MeshSampler<MyMesh> mps(subM);
    tri::MeshSampler<MyMesh> mrs(rndM);

    if (tri::io::Importer<MyMesh>::Open(m, argv[1]) != 0)
    {
        printf("Error reading file  %s\n", argv[1]);
        exit(0);
    }
    tri::SurfaceSampling<MyMesh, tri::TrivialSampler<MyMesh> >::SamplingRandomGenerator().initialize(time(0));
    float perc = atof(argv[2]);

    float radius = m.bbox.Diag() * perc;
    printf("Subsampling a PointCloud of %i vert with %f radius\n", m.VN(), radius);
    tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh> >::PoissonDiskParam pp;
    pp.bestSampleChoiceFlag = false;
    tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh> >::PoissonDiskPruning(mps, m, radius, pp);
    tri::io::ExporterPLY<MyMesh>::Save(subM, "PoissonMesh.ply");
    printf("Sampled %i vertices in %5.2f\n", subM.VN(), float(pp.pds.pruneTime + pp.pds.gridTime) / float(CLOCKS_PER_SEC));

    int t0 = clock();
    tri::Clustering<MyMesh, vcg::tri::AverageColorCell<MyMesh> > ClusteringGrid;
    ClusteringGrid.Init(m.bbox, 10000, radius * sqrt(2.0f));
    ClusteringGrid.AddPointSet(m);
    ClusteringGrid.ExtractMesh(cluM);
    int t1 = clock();
    tri::io::ExporterPLY<MyMesh>::Save(cluM, "ClusterMesh.ply");
    printf("Sampled %i vertices in %5.2f\n", cluM.VN(), float(t1 - t0) / float(CLOCKS_PER_SEC));

    int t2 = clock();
    int sampleNum = (cluM.VN() + subM.VN()) / 2;
    tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh> >::VertexUniform(m, mrs, sampleNum);
    int t3 = clock();
    tri::io::ExporterPLY<MyMesh>::Save(rndM, "RandomMesh.ply");
    printf("Sampled %i vertices in %5.2f\n", rndM.VN(), float(t3 - t2) / float(CLOCKS_PER_SEC));


    return 0;
}*/


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