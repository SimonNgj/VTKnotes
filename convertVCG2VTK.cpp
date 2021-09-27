#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

#include <iostream>
#include <pcl/io/ply_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>

#include<vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <vtkPLYWriter.h>
#include <vtkPLYReader.h>
#include <vtkSTLReader.h>
#include <pcl/visualization/cloud_viewer.h>
#include <vtkNamedColors.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkTriangle.h>
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkCellArrayIterator.h"
#include "vtkPolyData.h"

using namespace vcg;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
    Use<MyFace>			::AsFaceType> {};

class MyVertex : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags, vertex::Mark> {};
class MyFace : public Face  < MyUsedTypes, face::VertexRef, face::Normal3f, face::BitFlags > {};
class MyMesh : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> > {};

//// Convert a PCL point cloud to VTK
//void PCL_TO_VTK(pcl::PointCloud<pcl::PointXYZRGB>& pcl_PointCloud, vtkSmartPointer<vtkPoints> vtk_PointCloud)
//{
//    for (int i = 0; i < (int)pcl_PointCloud.size(); i++)
//    {
//        double point[3];
//        point[0] = pcl_PointCloud.points[i].x;
//        point[1] = pcl_PointCloud.points[i].y;
//        point[2] = pcl_PointCloud.points[i].z;
//        vtk_PointCloud->InsertNextPoint(point);
//    }
//}
//
//// Convert a PCL mesh to a VCG mesh
//void Mesh_PCL2VCG(pcl::PolygonMesh mesh_binary2, MyMesh *m)
//{
//    m->Clear();
//    // Now convert the vertices to VCG MyMesh
//    int vertCount = mesh_binary2.cloud.width * mesh_binary2.cloud.height;
//    vcg::tri::Allocator<MyMesh>::AddVertices(*m, vertCount);
//    pcl::PointCloud<pcl::PointXYZ> cloud1;
//    pcl::fromPCLPointCloud2(mesh_binary2.cloud, cloud1);
//    for (unsigned int i = 0; i < vertCount; ++i)
//        m->vert[i].P() = vcg::Point3f(cloud1.points[i].x, cloud1.points[i].y, cloud1.points[i].z);
//
//    // Now convert the polygon indices to VCG MyMesh => make VCG faces..
//    int triCount = mesh_binary2.polygons.size();
//    if (triCount == 1)
//    {
//        if (mesh_binary2.polygons[0].vertices[0] == 0 && mesh_binary2.polygons[0].vertices[1] == 0 && mesh_binary2.polygons[0].vertices[2] == 0)
//            triCount = 0;
//    }
//    Allocator<MyMesh>::AddFaces(*m, triCount);
//    for (unsigned int i = 0; i < triCount; ++i)
//    {
//        m->face[i].V(0) = &m->vert[mesh_binary2.polygons[i].vertices[0]];
//        m->face[i].V(1) = &m->vert[mesh_binary2.polygons[i].vertices[1]];
//        m->face[i].V(2) = &m->vert[mesh_binary2.polygons[i].vertices[2]];
//    }
//
//    vcg::tri::UpdateBounding<MyMesh>::Box(*m);
//    vcg::tri::UpdateNormal<MyMesh>::PerFace(*m);
//    vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(*m);
//}
//
//// Convert a VCG mesh to a PCL mesh
//void Mesh_VCG2PCL(MyMesh *m, pcl::PolygonMesh& mesh_binary2)
//{
//    int vertCount = m->vert.size();
//    int triCount = m->face.size();
//
//    pcl::PointCloud<pcl::PointXYZ> cloud2;
//    for (int i = 0; i < vertCount; i++)
//    {
//        pcl::PointXYZ point;
//        point.x = m->vert[i].P()[0];
//        point.y = m->vert[i].P()[1];
//        point.z = m->vert[i].P()[2];
//        cloud2.push_back(point);
//    }
//
//    std::vector<pcl::Vertices> polygons2;
//
//    // Now fill the indices of the triangles/faces of the mesh
//    for (int i = 0; i < triCount; i++)
//    {
//        pcl::Vertices vertices2;
//        vertices2.vertices.push_back(m->face[i].V(0) - &*m->vert.begin());
//        vertices2.vertices.push_back(m->face[i].V(1) - &*m->vert.begin());
//        vertices2.vertices.push_back(m->face[i].V(2) - &*m->vert.begin());
//        polygons2.push_back(vertices2);
//    }
//    mesh_binary2.polygons = polygons2;
//    pcl::toPCLPointCloud2(cloud2, mesh_binary2.cloud);
//}

int main()
{
    ///////// Convert VTK mesh to VCG mesh ////////
    vtkNew<vtkPLYReader> reader;
    reader->SetFileName("../../files/10.ply");
    reader->Update();
    vtkSmartPointer<vtkPolyData> vtk_mesh = vtkSmartPointer<vtkPolyData>::New();
    vtk_mesh = reader->GetOutput();

    vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
    points1 = vtk_mesh->GetPoints();
    int vertCount = points1->GetNumberOfPoints();
    vtkSmartPointer<vtkCellArray> polys1 = vtkSmartPointer<vtkCellArray>::New();
    polys1 = vtk_mesh->GetPolys();
    int triCount = polys1->GetNumberOfCells();

    MyMesh m;
    m.Clear();
    // Now convert the vertices to VCG MyMesh
    vcg::tri::Allocator<MyMesh>::AddVertices(m, vertCount);
    for (unsigned int i = 0; i < vertCount; ++i)
    {
        double p[3];
        points1->GetPoint(i, p);

        m.vert[i].P() = vcg::Point3f(p[0], p[1], p[2]);
    }

    Allocator<MyMesh>::AddFaces(m, triCount);
    vtkNew<vtkIdList> idL;
    vtk_mesh->GetPolys()->InitTraversal();
    int id_poly = 0;
    while (vtk_mesh->GetPolys()->GetNextCell(idL))
    {
        m.face[id_poly].V(0) = &m.vert[idL->GetId(0)];
        m.face[id_poly].V(1) = &m.vert[idL->GetId(1)];
        m.face[id_poly].V(2) = &m.vert[idL->GetId(2)];
        id_poly++;
    }

    vcg::tri::UpdateBounding<MyMesh>::Box(m);
    vcg::tri::UpdateNormal<MyMesh>::PerFace(m);
    vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(m);

    std::string meshFilename1 = "../../files/10_VTK2VCG.ply";
    vcg::tri::io::Exporter<MyMesh>::Save(m, meshFilename1.c_str());

    /////////// Convert VTK mesh to VCG mesh ///////////

    ///////// Convert VCG mesh to VTK mesh ////////
    MyMesh m2;
    bool openFile2 = tri::io::Importer<MyMesh>::Open(m2, "../../files/10.ply");

    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();

    // Convert points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i<m2.vert.size(); i++)
    {
        points->InsertNextPoint(m2.vert[i].P()[0] , m2.vert[i].P()[1], m2.vert[i].P()[2]);
    }
    output->SetPoints(points);

    // Convert polys
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 0; i<m2.face.size(); i++)
    {
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

        triangle->GetPointIds()->SetId(0, m2.face[i].V(0) - &*m2.vert.begin());
        triangle->GetPointIds()->SetId(1, m2.face[i].V(1) - &*m2.vert.begin());
        triangle->GetPointIds()->SetId(2, m2.face[i].V(2) - &*m2.vert.begin());

        polys->InsertNextCell(triangle);
    }
    output->SetPolys(polys);

    // Write to an *.ply file
    vtkNew<vtkPLYWriter> plyWriter;
    plyWriter->SetFileName("../../files/10_VCG2VTK.ply");
    plyWriter->SetInputData(output);
    plyWriter->Write();
    ///////// Convert VCG mesh to VTK mesh ////////



    //pcl::PolygonMesh mesh_bin;
    //pcl::io::loadPLYFile("../../files/10.ply", mesh_bin);

    ///////////// Convert PCL 2 VCG ///////////
    //MyMesh *me = new MyMesh;
    //Mesh_PCL2VCG(mesh_bin, me);

    //vcg::tri::io::Exporter<MyMesh>::Save(*me, "../../files/10_testPCL2VCG.ply");
    ///////////// Convert PCL 2 VCG ///////////

    ///////////// Convert VCG 2 PCL ///////////
    //MyMesh* me2 = new MyMesh;
    //bool openFile = tri::io::Importer<MyMesh>::Open(*me2, "../../files/10.ply");
    //pcl::PolygonMesh mesh_pcl2;
    //Mesh_VCG2PCL(me2, mesh_pcl2);

    //pcl::io::savePLYFile("../../files/10_testVCG2PCL2.ply", mesh_pcl2);
    ///////////// Convert VCG 2 PCL ///////////

    return EXIT_SUCCESS;
}

/*
int main()
{
    //pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
    //pcl::io::loadPCDFile("../../files/car6.pcd", *cloud);
    pcl::PolygonMesh mesh_binary;
    pcl::io::loadPLYFile("../../files/10.ply", mesh_binary);

    
    /////////// Convert PCL 2 VCG ///////////
    MyMesh m;
    m.Clear();
    // Now convert the vertices to VCG MyMesh
    int vertCount = mesh_binary.cloud.width * mesh_binary.cloud.height;
    vcg::tri::Allocator<MyMesh>::AddVertices(m, vertCount);
    pcl::PointCloud<pcl::PointXYZ> cloud1;
    pcl::fromPCLPointCloud2(mesh_binary.cloud, cloud1);
    for (unsigned int i = 0; i < vertCount; ++i)
        m.vert[i].P() = vcg::Point3f(cloud1.points[i].x, cloud1.points[i].y, cloud1.points[i].z);

    // Now convert the polygon indices to VCG MyMesh => make VCG faces..
    int triCount = mesh_binary.polygons.size();
    if (triCount == 1)
    {
        if (mesh_binary.polygons[0].vertices[0] == 0 && mesh_binary.polygons[0].vertices[1] == 0 && mesh_binary.polygons[0].vertices[2] == 0)
            triCount = 0;
    }
    Allocator<MyMesh>::AddFaces(m, triCount);
    for (unsigned int i = 0; i < triCount; ++i)
    {
        m.face[i].V(0) = &m.vert[mesh_binary.polygons[i].vertices[0]];
        m.face[i].V(1) = &m.vert[mesh_binary.polygons[i].vertices[1]];
        m.face[i].V(2) = &m.vert[mesh_binary.polygons[i].vertices[2]];
    }

    vcg::tri::UpdateBounding<MyMesh>::Box(m);
    vcg::tri::UpdateNormal<MyMesh>::PerFace(m);
    vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(m);

    std::string meshFilename2 = "../../files/10_testPCL2VCG.ply";
    vcg::tri::io::Exporter<MyMesh>::Save(m, meshFilename2.c_str());

    /////////// Convert PCL 2 VCG ///////////

    /////////// Convert VCG 2 PCL ///////////
    pcl::PolygonMesh mesh_binary2;

    pcl::PointCloud<pcl::PointXYZ> cloud2;
    for (int i = 0; i < vertCount; i++)
    {
        pcl::PointXYZ point;
        point.x = m.vert[i].P()[0];
        point.y = m.vert[i].P()[1];
        point.z = m.vert[i].P()[2];
        cloud2.push_back(point);
    }


    std::vector<pcl::Vertices> polygons2;

    // Now fill the indices of the triangles/faces of the mesh
    for (int i = 0; i < triCount; i++)
    {
        pcl::Vertices vertices2;
        vertices2.vertices.push_back(m.face[i].V(0) - &*m.vert.begin());
        vertices2.vertices.push_back(m.face[i].V(1) - &*m.vert.begin());
        vertices2.vertices.push_back(m.face[i].V(2) - &*m.vert.begin());
        polygons2.push_back(vertices2);
    }
    mesh_binary2.polygons = polygons2;
    pcl::toPCLPointCloud2(cloud2, mesh_binary2.cloud);

    pcl::io::savePLYFile("../../files/10_testVCG2PCL.ply", mesh_binary2);

    /////////// Convert VCG 2 PCL ///////////


    /////////// Visualization ///////////
    //vtkSmartPointer<vtkPoints> output = vtkSmartPointer<vtkPoints>::New();
    //PCL_TO_VTK(*cloud, output);

    //vtkNew<vtkPolyData> polyData;

    //polyData->SetPoints(output);
    /////////// Convert PCL 2 VTK///////////

    //vtkNew<vtkVertexGlyphFilter> glyphFilter;
    //glyphFilter->SetInputData(polyData);
    //glyphFilter->Update();

    //// Visualize
    //vtkNew<vtkNamedColors> colors;
    //vtkNew<vtkPolyDataMapper> mapper;
    //mapper->SetInputConnection(glyphFilter->GetOutputPort());

    //vtkNew<vtkActor> actor;
    //actor->SetMapper(mapper);
    //actor->GetProperty()->SetPointSize(2);
    //actor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

    //vtkNew<vtkRenderer> renderer;
    //renderer->AddActor(actor);
    //renderer->SetBackground(colors->GetColor3d("Gainsboro").GetData());

    //vtkNew<vtkRenderWindow> renderWindow;
    //renderWindow->AddRenderer(renderer);
    //renderWindow->SetWindowName("ReadTextFile");

    //vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    //renderWindowInteractor->SetRenderWindow(renderWindow);
    //vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    //style->SetCurrentRenderer(renderer);
    //renderWindow->GetInteractor()->SetInteractorStyle(style);

    //renderWindow->Render();
    //renderWindowInteractor->Start();

    return EXIT_SUCCESS;
    //////////////////////
}*/