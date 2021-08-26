//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2)
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//
//#include "subFunc.h"
//
//int main(int argc, char** argv)
//{
//
//    MyMesh m;
//
//    std::string namefile = "../files/arrow";
//    bool openFile = tri::io::Importer<MyMesh>::Open(m, (namefile + ".off").c_str());
//    //bool openFile = tri::io::Importer<MyMesh>::Open(m, argv[1]);
//
//    ////--------------------  Calculate the normals  -------------------------
//    tri::PointCloudNormal<MyMesh>::Param p;
//    p.fittingAdjNum = 10;
//    p.smoothingIterNum = 0;
//    tri::PointCloudNormal<MyMesh>::Compute(m, p, 0);
//
//    //vcg::tri::io::Exporter<MyMesh>::Save(m, "../files/10_calNormal.off");
//    ////------------------ End calculate the normals ------------------------
//
//    ////-------------------- Screened Poisson -----------------------------
//    // Convert to PCL format
//    PointCloud<PointNormal>::Ptr cloud_smoothed_normals(new PointCloud<PointNormal>());
//
//    for (auto vp = m.vert.begin(); vp != m.vert.end(); ++vp)
//    {
//        if (!vp->IsD())
//        {
//            pcl::PointNormal point;
//            point.x = vp->P()[0];
//            point.y = vp->P()[1];
//            point.z = vp->P()[2];
//            point.normal_x = vp->N()[0];
//            point.normal_y = vp->N()[1];
//            point.normal_z = vp->N()[2];
//            cloud_smoothed_normals->push_back(point);
//        }
//    }
//
//    // Apply screened poisson function
//    cout << "Begin poisson reconstruction..." << endl;
//    Poisson<PointNormal> poisson;
//    poisson.setDepth(8);
//    poisson.setPointWeight(4);
//    poisson.setScale(1.1);
//    poisson.setInputCloud(cloud_smoothed_normals);
//    PolygonMesh mesh;
//    poisson.reconstruct(mesh);
//
//    //string meshFilename = "..\\files\\10_poisson.ply";
//    //pcl::io::savePLYFile(meshFilename, mesh);
//    string meshFilename = namefile + "_poisson.stl";
//    pcl::io::savePolygonFileSTL(meshFilename, mesh);
//
//    //////////////////// Read STL back //////////////////////////
//    MyMesh m2;
//    //MyMesh m3;
//    
//    bool openFile2 = tri::io::Importer<MyMesh>::Open(m2, meshFilename.c_str());
//
//    ///////////// Remove duplicate vertices ///////////////////
//    int delvert = tri::Clean<MyMesh>::RemoveDuplicateVertex(m2);
//    tri::Clean<MyMesh>::RemoveUnreferencedVertex(m2);
//    Allocator<MyMesh>::CompactEveryVector(m2);
//    //tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(m2);
//    //tri::UpdateBounding<MyMesh>::Box(m2);
//
//    //////////////// Merge close vertices /////////////////////
//    float threshold = 2.0;
//    int total = tri::Clean<MyMesh>::MergeCloseVertex(m2, threshold);
//    tri::Clean<MyMesh>::RemoveUnreferencedVertex(m2);
//    Allocator<MyMesh>::CompactEveryVector(m2);
//
//    //////////////// Remove Isolated pieces (wrt Face Num.) /////////////////////
//    int minCC = 25;
//    std::pair<int, int> delInfo = tri::Clean<MyMesh>::RemoveSmallConnectedComponentsSize(m2, minCC);
//
//    //"removeUnref"
//    int delvert = tri::Clean<MyMesh>::RemoveUnreferencedVertex(m2);
//    tri::Clean<MyMesh>::RemoveUnreferencedVertex(m2);
//    Allocator<MyMesh>::CompactEveryVector(m2);
//
//    /// ////////////// Save /////////////////////
//    string meshFilename2 = namefile + "_poisson_refine.stl";
//    vcg::tri::io::Exporter<MyMesh>::Save(m2, meshFilename2.c_str());
//    //vcg::tri::io::Exporter<MyMesh>::Save(m2, "../files/remesh.ply");
//
//    //////--------------------- End Screened Poisson ----------------------------
//    //////////////////
//    vtkSmartPointer<vtkPolyData> vtk_mesh = vtkSmartPointer<vtkPolyData>::New();
//    pcl::io::mesh2vtk(mesh, vtk_mesh);
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(vtk_mesh);
//    triangleFilter->Update();
//
//    vtkNew<vtkPolyDataMapper> mapper2;
//    mapper2->SetInputConnection(triangleFilter->GetOutputPort());
//
//    //////////////////
//    //vtkNew<vtkPLYReader> reader;
//    //reader->SetFileName(meshFilename.c_str());
//    //vtkNew<vtkPolyDataMapper> mapper2;
//    //mapper2->SetInputConnection(reader->GetOutputPort());
//
//    vtkNew<vtkNamedColors> colors;
//    vtkNew<vtkActor> actor2;
//    actor2->SetMapper(mapper2);
//    //actor2->GetProperty()->SetRepresentationToWireframe();
//    //actor2->GetProperty()->ShadingOff();
//    actor2->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
//
//    /////////////
//    // Convert a VCG point cloud to vtk
//    vtkSmartPointer<vtkPoints> vtk_input = vtkSmartPointer<vtkPoints>::New();
//
//    for (auto vp = m.vert.begin(); vp != m.vert.end(); ++vp)
//    {
//        if (!vp->IsD())
//        {
//            double point[3];
//            point[0] = vp->P()[0];
//            point[1] = vp->P()[1];
//            point[2] = vp->P()[2];
//            vtk_input->InsertNextPoint(point);
//        }
//    }
//
//    vtkNew<vtkPolyData> polyData;
//
//    polyData->SetPoints(vtk_input);
//
//    vtkNew<vtkVertexGlyphFilter> glyphFilter;
//    glyphFilter->SetInputData(polyData);
//    glyphFilter->Update();
//
//    // Visualize
//    vtkNew<vtkPolyDataMapper> mapper1;
//    mapper1->SetInputConnection(glyphFilter->GetOutputPort());
//
//
//    vtkNew<vtkActor> actor1;
//    actor1->SetMapper(mapper1);
//    actor1->GetProperty()->SetPointSize(2);
//    actor1->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
//
//    vtkNew<vtkRenderer> renderer1;
//    renderer1->SetViewport(0., 0., 0.5, 1.);
//    //renderer1->AddActor(actor);
//    renderer1->SetBackground(colors->GetColor3d("Gainsboro").GetData());
//
//    vtkNew<vtkRenderer> renderer2;
//    renderer2->SetViewport(0.5, 0., 1., 1.);
//    renderer2->SetBackground(colors->GetColor3d("Gainsboro").GetData());
//
//    vtkNew<vtkRenderWindow> renderWindow;
//    renderWindow->SetSize(800, 400);
//    renderWindow->AddRenderer(renderer1);
//    renderWindow->AddRenderer(renderer2);
//    renderWindow->SetWindowName("PointCloud2Mesh");
//
//    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
//    style->SetCurrentRenderer(renderer1);
//    renderWindow->GetInteractor()->SetInteractorStyle(style);
//
//    renderer1->AddViewProp(actor1);
//    renderer2->AddViewProp(actor2);
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    //system("pause");
//
//    return 0;
//}