#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
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

int main()
{
    pcl::PolygonMesh mesh_bin;
    pcl::io::loadPLYFile("../../files/10.ply", mesh_bin);

    pcl::PointCloud<pcl::PointXYZ> cloud1;
    pcl::fromPCLPointCloud2(mesh_bin.cloud, cloud1);
    int vertCount = cloud1.width * cloud1.height;
    int triCount = mesh_bin.polygons.size();

    ///////// Convert PCL mesh to VTK mesh ////////
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();

    // Convert points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto pt : cloud1.points) 
    {
        points->InsertNextPoint(pt.x, pt.y, pt.z);
    }
    output->SetPoints(points);

    // Convert polys
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    for (auto t : mesh_bin.polygons) 
    {
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

        triangle->GetPointIds()->SetId(0, t.vertices[0]);
        triangle->GetPointIds()->SetId(1, t.vertices[1]);
        triangle->GetPointIds()->SetId(2, t.vertices[2]);

        polys->InsertNextCell(triangle);
    }
    output->SetPolys(polys);

    // Write to an *.ply file
    vtkNew<vtkPLYWriter> plyWriter;
    plyWriter->SetFileName("../../files/10_PCL2VTK.ply");
    plyWriter->SetInputData(output);
    plyWriter->Write();

    ///////// Convert VTK mesh to PCL mesh ////////
    vtkNew<vtkPLYReader> reader;
    reader->SetFileName("../../files/10_PCL2VTK.ply");
    reader->Update();
    vtkSmartPointer<vtkPolyData> vtk_mesh = vtkSmartPointer<vtkPolyData>::New();
    vtk_mesh = reader->GetOutput();

    pcl::PolygonMesh mesh_binary2;

    vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
    points2 = vtk_mesh->GetPoints();
    int noV = points2->GetNumberOfPoints();
    vtkSmartPointer<vtkCellArray> polys2 = vtkSmartPointer<vtkCellArray>::New();
    polys2 = vtk_mesh->GetPolys();
    int noP = polys2->GetNumberOfCells();
   
    pcl::PointCloud<pcl::PointXYZ> cloud2;
    for (int i = 0; i < points2->GetNumberOfPoints(); i++)
    {
        pcl::PointXYZ point;
        double p[3];
        points2->GetPoint(i, p);
        
        point.x = p[0];
        point.y = p[1];
        point.z = p[2];
        cloud2.push_back(point);
    }
    pcl::toPCLPointCloud2(cloud2, mesh_binary2.cloud);

    std::vector<pcl::Vertices> polygons2;

    //int id_poly = 0;
    vtkNew<vtkIdList> idL;
    vtk_mesh->GetPolys()->InitTraversal();
    while (vtk_mesh->GetPolys()->GetNextCell(idL))
    {
        pcl::Vertices vertices2;

        vertices2.vertices.push_back(idL->GetId(0));
        vertices2.vertices.push_back(idL->GetId(1));
        vertices2.vertices.push_back(idL->GetId(2));

        polygons2.push_back(vertices2);
    }
    mesh_binary2.polygons = polygons2;

    pcl::io::savePLYFile("../../files/10_VTK2PCL.ply", mesh_binary2);

   /* pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::io::loadPCDFile("..\\..\\files\\car6.pcd", *cloud);*/

    ////////////////////////
    //vtkSmartPointer<vtkPoints> output = vtkSmartPointer<vtkPoints>::New();
    //PCL_TO_VTK(*cloud, output);

 /*   vtkNew<vtkPolyData> polyData;

    polyData->SetPoints(output);*/

    //vtkNew<vtkVertexGlyphFilter> glyphFilter;
    //glyphFilter->SetInputData(vtk_mesh);
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
}