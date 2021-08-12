#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

#include <vtkActor.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkNamedColors.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtksys/SystemTools.hxx>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkTransform.h>
#include <vtkPLYReader.h>
#include <vtkTransformPolyDataFilter.h>

#include "subFunc.h"

vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName)
{
	vtkSmartPointer<vtkPolyData> polyData;
	std::string extension = vtksys::SystemTools::GetFilenameExtension(std::string(fileName));

	if (extension == ".STL")
	{
		vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(fileName);
		reader->Update();
		polyData = reader->GetOutput();
	}
	else
	{
		vtkSmartPointer<vtkSphereSource> source = vtkSmartPointer<vtkSphereSource>::New();
		source->Update();
		polyData = source->GetOutput();
	}
	return polyData;
}

int main(int, char *[])
{
	vtkSmartPointer<vtkPolyData> poly1;
	vtkSmartPointer<vtkPolyData> poly2;

	string in1 = "../files/exp2buildobj.STL";
	string in2 = "../files/exp2sub1smooth.STL";
	poly1 = ReadPolyData(in1.c_str());
	poly2 = ReadPolyData(in2.c_str());

	vtkSmartPointer<vtkTransform> translation1 = vtkSmartPointer<vtkTransform>::New();
	translation1->Translate(60, 20, -30);
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter1 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter1->SetInputData(poly1);
	transformFilter1->SetTransform(translation1);
	transformFilter1->Update();

	//vtkSmartPointer<vtkTransform> translation2 = vtkSmartPointer<vtkTransform>::New();
	//translation2->RotateX(180);
	//vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter2 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	//transformFilter2->SetInputData(poly2);
	//transformFilter2->SetTransform(translation2);
	//transformFilter2->Update();

	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

	vtkSmartPointer<vtkPolyDataMapper> input1Mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	input1Mapper->SetInputData(transformFilter1->GetOutput());
	vtkSmartPointer<vtkActor> input1Actor = vtkSmartPointer<vtkActor>::New();
	input1Actor->GetProperty()->SetOpacity(0.8);
	input1Actor->SetMapper(input1Mapper);
	input1Actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Red").GetData());
	//input1Actor->SetPosition(0, 0, 0);

	vtkSmartPointer<vtkPolyDataMapper> input2Mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	input2Mapper->SetInputData(poly2);
	vtkSmartPointer<vtkActor> input2Actor = vtkSmartPointer<vtkActor>::New();
	//input2Actor->GetProperty()->SetOpacity(0.2);
	input2Actor->SetMapper(input2Mapper);
	input2Actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Blue").GetData());
	//input2Actor->SetPosition(0, 0, 0);

	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
	//booleanOperation->SetOperationToUnion();
	//booleanOperation->SetOperationToIntersection();
	booleanOperation->SetOperationToDifference();
	booleanOperation->SetInputData(0, transformFilter1->GetOutput());
	booleanOperation->SetInputData(1, poly2); // poly2

	vtkNew<vtkSTLWriter> stlWriter;
	string out1 = in1 + "compare.STL";
	stlWriter->SetFileName(out1.c_str());
	stlWriter->SetInputConnection(booleanOperation->GetOutputPort());
	stlWriter->Write();

	vtkSmartPointer<vtkPolyDataMapper> booleanOperationMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	booleanOperationMapper->SetInputConnection(booleanOperation->GetOutputPort());
	booleanOperationMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> booleanOperationActor = vtkSmartPointer<vtkActor>::New();
	//booleanOperationActor->GetProperty()->SetOpacity(0.5);
	booleanOperationActor->SetMapper(booleanOperationMapper);
	booleanOperationActor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Green").GetData());
	booleanOperationActor->GetProperty()->SetSpecular(.6);
	booleanOperationActor->GetProperty()->SetSpecularPower(20);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	//renderer->AddViewProp(input1Actor);
	//renderer->AddViewProp(input2Actor);
	renderer->AddViewProp(booleanOperationActor);

	renderer->SetBackground(colors->GetColor3d("Black").GetData());
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(640, 480);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renderWindow);
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    style->SetCurrentRenderer(renderer);
    renderWindow->GetInteractor()->SetInteractorStyle(style);

	renderWindow->Render();
	renWinInteractor->Start();

	return EXIT_SUCCESS;
}
