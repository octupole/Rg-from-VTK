/*
 * StructuredGrid.cpp
 *
 *  Created on: Dec 3, 2018
 *      Author: marchi
 */
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include "vtkFloatArray.h"
#include "vtkDataArray.h"
#include "vtkAbstractArray.h"

#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include "vtkPointData.h"

#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

const int XX{0},YY{1},ZZ{2},DIM{3};

int main(int, char *[])
{
  // Create a grid
	std::string inputFilename ="ShuA-charmm.vts";
	vtkSmartPointer<vtkStructuredGrid> sGrid =vtkSmartPointer<vtkStructuredGrid>::New();
	vtkSmartPointer<vtkFloatArray> Density=vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkXMLStructuredGridReader> reader =vtkSmartPointer<vtkXMLStructuredGridReader>::New();

	reader->SetFileName(inputFilename.c_str());
	reader->Update();

	sGrid=reader->GetOutput();
	int dims[DIM];
	sGrid->GetDimensions(dims);
	double point[DIM];
	size_t tot{0},N{0};

	float c0[DIM],norm{0};


	int Size=dims[XX]*dims[YY]*dims[ZZ];
	Density->SetNumberOfTuples(Size);
	Density->InsertTuples(0,Size,0,sGrid->GetPointData()->GetAbstractArray(0));


	for(size_t o{0};o<dims[0];o++)
		for(size_t p{0};p<dims[0];p++)
			for(size_t q{0};q<dims[0];q++){
				sGrid->GetPoint(o,p,q,point);
				float density=Density->GetTuple1(N);
				c0[XX]+=point[XX]*density;
				c0[YY]+=point[YY]*density;
				c0[ZZ]+=point[ZZ]*density;
				norm+=density;
				N++;
			}
	for(size_t o{0};o<3;o++)
		c0[o]/=norm;
	float Rg;
	float x0,y0,z0;
	N=0;
	for(size_t o{0};o<dims[0];o++)
		for(size_t p{0};p<dims[0];p++)
			for(size_t q{0};q<dims[0];q++){
				sGrid->GetPoint(o,p,q,point);
				float density=Density->GetTuple1(N);
				x0=point[XX]-c0[XX];
				y0=point[YY]-c0[YY];
				z0=point[ZZ]-c0[ZZ];

				Rg+=(x0*x0+y0*y0+z0*z0)*density;
				N++;
			}
		Rg/=norm;





	//cout << lula->GetArrayType() << endl;

	cout << sqrt(Rg) << endl;

	//structuredGrid->GetPoint(1,2,3)

//	vtkSmartPointer<vtkStructuredGridReader> reader=vtkSmartPointer<vtkStructuredGridReader>::New();
//	reader->SetFileName(inputFilename.c_str());
//	reader->Update();

//  double x, y, z;
//
//  x = 0.0;
//  y = 0.0;
//  z = 0.0;
//
//  for(unsigned int k = 0; k < 2; k++)
//    {
//    z += 2.0;
//    for(unsigned int j = 0; j < 3; j++)
//      {
//      y += 1.0;
//      for(unsigned int i = 0; i < 2; i++)
//        {
//        x += .5;
//        points->InsertNextPoint(x, y, z);
//        }
//      }
//    }
//
//  // Specify the dimensions of the grid
//  structuredGrid->SetDimensions(2,3,2);
//  structuredGrid->SetPoints(points);
//
//  int* dims = structuredGrid->GetDimensions();
//
//  // Retrieve the entries from the grid and print them to the screen
//  unsigned int counter = 0;
//
//  for (int k = 0; k < dims[2]; k++)
//    {
//    for (int j = 0; j < dims[1]; j++)
//      {
//      for (int i = 0; i < dims[0]; i++)
//        {
//        double p[3];
//        structuredGrid->GetPoint(counter, p);
//
//        double pNew[3];
//        structuredGrid->GetPoint(i, j, k, pNew);
//
//        std::cout << "P   : "
//                  << p[0] << " "
//                  << p[1] << " "
//                  << p[2] << std::endl;
//        std::cout << "PNew: "
//                  << pNew[0] << " "
//                  << pNew[1] << " "
//                  << pNew[2] << std::endl;
//
//        counter++;
//        }
//      }
//    }
//
//  // Create a mapper and actor
//  vtkSmartPointer<vtkDataSetMapper> mapper =
//    vtkSmartPointer<vtkDataSetMapper>::New();
//  mapper->SetInputData(structuredGrid);
//
//  vtkSmartPointer<vtkActor> actor =
//    vtkSmartPointer<vtkActor>::New();
//  actor->SetMapper(mapper);
//
//  // Create a renderer, render window, and interactor
//  vtkSmartPointer<vtkRenderer> renderer =
//    vtkSmartPointer<vtkRenderer>::New();
//  vtkSmartPointer<vtkRenderWindow> renderWindow =
//    vtkSmartPointer<vtkRenderWindow>::New();
//  renderWindow->AddRenderer(renderer);
//  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//    vtkSmartPointer<vtkRenderWindowInteractor>::New();
//  renderWindowInteractor->SetRenderWindow(renderWindow);
//
//  // Add the actor to the scene
//  renderer->AddActor(actor);
//  renderer->SetBackground(.3, .6, .3); // Background color green
//
//  // Render and interact
//  renderWindow->Render();
//  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}




