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
#include <sstream>

const int XX{0},YY{1},ZZ{2},DIM{3};

int main(int argc, char *argv[])
{
	ifstream ftest;

	std::string inputFilename ="";
	try{
		std::stringstream ss;
		ss<<argv[0];
		if(argc == 3){
			if(strcmp (argv[1],"-i") == 0){
				inputFilename=argv[2];
				ftest.open(inputFilename.c_str(),ios::in);
				if(!ftest) throw std::string("\n Cannot open file: " + inputFilename + "\n");
				ftest.close();
			} else{
				throw std::string("Usage:\t"+ss.str() + " -i input .vts file");
			}
		}else {
			throw std::string("Usage:\t"+ss.str() + " -i input .vts file");
		}

	}
	catch(const std::string & s){
		cout << s << endl;
		exit(1);
	}

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
	cout << sqrt(Rg) << endl;
  return EXIT_SUCCESS;
}




