#include "Grid.h"
#include "PathwaySet.h"

void Grid::write_K_field_to_vtp(string filename, double z_factor)
{
	vtkSmartPointer<vtkPoints> points_3 =
		vtkSmartPointer<vtkPoints>::New();

	double xx, yy, zz;
	vtkSmartPointer<vtkFloatArray> values =
		vtkSmartPointer<vtkFloatArray>::New();

	values->SetNumberOfComponents(1);

	values->SetName("K");


	for (unsigned int x = 0; x < geometric_parameters.nx; x++)
	{
		for (unsigned int y = 0; y < geometric_parameters.ny; y++)
		{
			xx = x*geometric_parameters.dx+unitrandom()*0.05*geometric_parameters.dx;
			yy = y*geometric_parameters.dy+unitrandom()*0.05*geometric_parameters.dy;
			zz = cell_property[x][y].K*z_factor;


            float t[1] = { float(cell_property[x][y].K) };
            points_3->InsertNextPoint(xx, yy, zz);
            values->InsertNextTupleValue(t);

		}
	}

	// Add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> inputPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	inputPolyData->SetPoints(points_3);

	// Triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay =
		vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
	delaunay->SetInput(inputPolyData);
#else
	delaunay->SetInputData(inputPolyData);
#endif
	delaunay->Update();
	vtkPolyData* outputPolyData = delaunay->GetOutput();

	double bounds[6];
	outputPolyData->GetBounds(bounds);

	// Find min and max z
	double minz = bounds[4];
	double maxz = bounds[5];

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		outputPolyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);
		//std::cout << "dcolor: "
		//	<< dcolor[0] << " "
		//	<< dcolor[1] << " "
		//	<< dcolor[2] << std::endl;
		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
		//std::cout << "color: "
		//	<< (int)color[0] << " "
		//	<< (int)color[1] << " "
		//	<< (int)color[2] << std::endl;

		colors_2->InsertNextTupleValue(color);
	}

	outputPolyData->GetPointData()->SetScalars(values);


	//Append the two meshes
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	//appendFilter->AddInputData(polydata);
	//appendFilter->AddInputData(polydata_1);
	appendFilter->AddInputData(outputPolyData);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();


}

void Grid::write_field_to_vtp(const string &quantity, const string &filename, double z_factor)
{
	vtkSmartPointer<vtkPoints> points_3 =
		vtkSmartPointer<vtkPoints>::New();

	double xx, yy, zz;
	vtkSmartPointer<vtkFloatArray> values =
		vtkSmartPointer<vtkFloatArray>::New();

	values->SetNumberOfComponents(1);

	values->SetName(quantity.c_str());


	for (unsigned int x = 0; x < geometric_parameters.nx; x++)
	{
		for (unsigned int y = 0; y < geometric_parameters.ny; y++)
		{
			xx = x*geometric_parameters.dx+unitrandom()*0.05*geometric_parameters.dx;
			yy = y*geometric_parameters.dy+unitrandom()*0.05*geometric_parameters.dy;
			zz = cell_property[x][y].PropertyValue(quantity)*z_factor;


            float t[1] = { float(cell_property[x][y].PropertyValue(quantity)) };
            points_3->InsertNextPoint(xx, yy, zz);
            values->InsertNextTupleValue(t);

		}
	}

	// Add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> inputPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	inputPolyData->SetPoints(points_3);

	// Triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay =
		vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
	delaunay->SetInput(inputPolyData);
#else
	delaunay->SetInputData(inputPolyData);
#endif
	delaunay->Update();
	vtkPolyData* outputPolyData = delaunay->GetOutput();

	double bounds[6];
	outputPolyData->GetBounds(bounds);

	// Find min and max z
	double minz = bounds[4];
	double maxz = bounds[5];

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		outputPolyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);
		//std::cout << "dcolor: "
		//	<< dcolor[0] << " "
		//	<< dcolor[1] << " "
		//	<< dcolor[2] << std::endl;
		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
		//std::cout << "color: "
		//	<< (int)color[0] << " "
		//	<< (int)color[1] << " "
		//	<< (int)color[2] << std::endl;

		colors_2->InsertNextTupleValue(color);
	}

	outputPolyData->GetPointData()->SetScalars(values);


	//Append the two meshes
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	//appendFilter->AddInputData(polydata);
	//appendFilter->AddInputData(polydata_1);
	appendFilter->AddInputData(outputPolyData);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();


}


void PathwaySet::SaveToVTP(vtkSmartPointer<vtkPolyDataMapper> mapper, string filename)
{


	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();

}

vtkSmartPointer<vtkPolyDataMapper> PathwaySet::ConvertToVTP(double z_factor, double offset)
{
	vector<vtkSmartPointer<vtkPolyData>> outarray;
	for (int i = 0; i < paths.size(); i++)
		outarray.push_back(paths[i].ToVTP(z_factor, offset, i));

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	for (int i = 0; i < outarray.size(); i++)
		appendFilter->AddInputData(outarray[i]);
#endif
	appendFilter->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
#endif

	return mapper;
}

vtkSmartPointer<vtkPolyData> Pathway::ToVTP(double z_factor, double offset, int ID)
{

	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> positions =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkFloatArray> values_ID =
		vtkSmartPointer<vtkFloatArray>::New();
	values_ID->SetNumberOfComponents(1);
	values_ID->SetName("ID");

	for (int i = 0; i<int(points.size()); i++)
	{
		double p[3] = { points[i].X() , points[i].Y(), offset };
		double _ID[1] = { float(ID) };
		positions->InsertNextPoint(p);
		values_ID->InsertNextTuple(_ID);

	}
	vtkSmartPointer<vtkPolyLine> polyLine =
		vtkSmartPointer<vtkPolyLine>::New();

	polyLine->GetPointIds()->SetNumberOfIds(points.size());
	for (unsigned int i = 0; i < points.size(); i++)
	{
		polyLine->GetPointIds()->SetId(i, i);
	}

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();

	// Add the points to the dataset
	polyData->SetPoints(positions);

	// Add the lines to the dataset
	polyData->SetLines(cells);

	polyData->GetPointData()->SetScalars(values_ID);
	//polyData->GetPointData()->AddArray(values_u);
	//polyData->GetPointData()->AddArray(values_z);
	//polyData->GetPointData()->AddArray(values_t);
	//polyData->GetPointData()->AddArray(values_v_eff);
	//polyData->GetPointData()->AddArray(values_t_eff);

	// Visualization

	return polyData;
}

void PathwaySet::SaveToVTP(string filename, double z_factor, double offset, int interval)
{
	vector<vtkSmartPointer<vtkPolyData>> outputmappers;

	//set_progress_value(0);
	for (int i = 0; i < paths.size(); i+=interval)
	{	outputmappers.push_back(paths[i].ToVTP(z_factor,offset,i));
        //set_progress_value((double)i/(double)Traj.n());
	}

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	for (int i=0; i<outputmappers.size(); i++)
		appendFilter->AddInputData(outputmappers[i]);
#endif
	appendFilter->Update();



	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	#ifdef QT_version
	main_window->get_ui()->ShowOutput->append("Writing vtp file... ");
	#endif // QT_version
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();

}
