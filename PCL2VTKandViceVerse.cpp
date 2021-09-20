#include <iostream>
#include <pcl/common/io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PolygonMesh.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>

using namespace std;

int main()
{
    //Load a pcl mesh
    pcl::PolygonMesh mesh;
    pcl::io::loadPolygonFilePLY("....ply", mesh);


    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());

    // Convert pcl mesh to vtk mesh
    pcl::io::mesh2vtk(mesh, polydata);
	
	// Convert a vtk mesh to pcl point cloud
    pcl::io::vtkPolyDataToPointCloud(polydata, *cloud);
    pcl::io::savePCDFileASCII("....pcd", *cloud);
    return 0;
}