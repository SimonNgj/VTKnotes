//#include "vtkAutoInit.h" 
//VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//
//#include <vtkActor.h>
//#include <vtkBox.h>
//#include <vtkCamera.h>
//#include <vtkColor.h>
//#include <vtkContourFilter.h>
//#include <vtkImplicitBoolean.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSampleFunction.h>
//#include <vtkSphere.h>
//#include <vtkTriangleFilter.h>
//#include <vtkCleanPolyData.h>
//#include <vtkPolyDataReader.h>
//#include <vtkSTLReader.h>
//#include <vtksys/SystemTools.hxx>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkSTLWriter.h>
//#include <vtkSphereSource.h>
//
//vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName)
//{
//	vtkSmartPointer<vtkPolyData> polyData;
//	std::string extension = vtksys::SystemTools::GetFilenameExtension(std::string(fileName));
//
//	if (extension == ".STL")
//	{
//		vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//		reader->SetFileName(fileName);
//		reader->Update();
//		polyData = reader->GetOutput();
//	}
//	else
//	{
//		vtkSmartPointer<vtkSphereSource> source = vtkSmartPointer<vtkSphereSource>::New();
//		source->Update();
//		polyData = source->GetOutput();
//	}
//	return polyData;
//}
//
//int main(int argc, char* argv[])
//{
//	// Define colors
//	vtkNew<vtkNamedColors> colors;
//	vtkColor3d actorColor = colors->GetColor3d("AliceBlue");
//	vtkColor3d EdgeColour = colors->GetColor3d("SteelBlue");
//	vtkColor3d BackgroundColour = colors->GetColor3d("Silver");
//
//	// create a sphere
//	vtkNew<vtkSphere> sphere;
//	sphere->SetCenter(1.0, 0.0, 0.0);
//	sphere->SetRadius(1);
//
//	// create a box
//	vtkNew<vtkBox> box;
//	box->SetBounds(-1, 1, -1, 1, -1, 1);
//
//
//	////vtkNew<vtkSTLReader> reader1;
//	////reader1->SetFileName("../files/block.STL");
//	////reader1->Update();
//
//	////vtkNew<vtkSTLReader> reader2;
//	////reader2->SetFileName("../files/arrow.STL");
//	////reader2->Update();
//
//	//vtkSmartPointer<vtkPolyData> input1;
//	//vtkSmartPointer<vtkPolyData> input2;
//
//	//vtkSmartPointer<vtkPolyData> poly1;
//	//vtkSmartPointer<vtkPolyData> poly2;
//
//	//poly1 = ReadPolyData("../files/block.STL");
//	//poly2 = ReadPolyData("../files/arrow.STL");
//
//	//vtkSmartPointer<vtkTriangleFilter> tri1 = vtkSmartPointer<vtkTriangleFilter>::New();
//	//tri1->SetInputData(poly1);
//	//vtkSmartPointer<vtkCleanPolyData> clean1 = vtkSmartPointer<vtkCleanPolyData>::New();
//	//clean1->SetInputConnection(tri1->GetOutputPort());
//	//clean1->Update();
//	//input1 = clean1->GetOutput();
//
//	//vtkSmartPointer<vtkTriangleFilter> tri2 = vtkSmartPointer<vtkTriangleFilter>::New();
//	//tri2->SetInputData(poly2);
//	//tri2->Update();
//	//vtkSmartPointer<vtkCleanPolyData> clean2 = vtkSmartPointer<vtkCleanPolyData>::New();
//	//clean2->SetInputConnection(tri2->GetOutputPort());
//	//clean2->Update();
//	//input2 = clean2->GetOutput();
//
//
//
//
//
//	// combine the two implicit functions
//	vtkSmartPointer<vtkImplicitBoolean> boolean;
//	boolean->SetOperationTypeToDifference();
//	//boolean->SetOperationTypeToUnion();
//	// boolean->SetOperationTypeToIntersection();
//	
//	//boolean->AddFunction(sphere);
//	//boolean->AddFunction(box);
//
//	boolean->AddFunction(box);
//	boolean->AddFunction(sphere);
//	
//	//boolean->AddFunction(reader2->GetOutputPort());
//
//
//	// The sample function generates a distance function from the implicit
//	// function.This is then contoured to get a polygonal surface.
//	vtkNew<vtkSampleFunction> sample;
//	sample->SetImplicitFunction(boolean);
//	sample->SetModelBounds(-1, 2, -1, 1, -1, 1);
//	sample->SetSampleDimensions(40, 40, 40);
//	sample->ComputeNormalsOff();
//
//	// contour
//	vtkNew<vtkContourFilter> surface;
//	surface->SetInputConnection(sample->GetOutputPort());
//	surface->SetValue(0, 0.0);
//
//	vtkNew<vtkSTLWriter> stlWriter;
//	stlWriter->SetFileName("../files/compareSTL.stl");
//	stlWriter->SetInputConnection(surface->GetOutputPort());
//	stlWriter->Write();
//
//	//////////////////////////
//
//	//////////////////////////
//
//	// Create a mapper and an actor
//	vtkNew<vtkPolyDataMapper> mapper;
//	mapper->SetInputConnection(surface->GetOutputPort());
//	mapper->ScalarVisibilityOff();
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->GetProperty()->EdgeVisibilityOn();
//	actor->GetProperty()->SetColor(actorColor.GetData());
//	actor->GetProperty()->SetEdgeColor(EdgeColour.GetData());
//
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	renderer->SetBackground(BackgroundColour.GetData());
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("BooleanOperationImplicitFunctions");
//
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// add the actor
//	renderer->AddActor(actor);
//
//	// Start
//	//renderer->GetActiveCamera()->SetPosition(5.0, -4.0, 1.6);
//	//renderer->GetActiveCamera()->SetViewUp(0.1, 0.5, 0.9);
//	//renderer->GetActiveCamera()->SetDistance(6.7);
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}