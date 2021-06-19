#include "Grid.h"
#include "NormalDist.h"
#define CMatrix_arma CMatrix


Grid::Grid()
{
    //ctor
}

Grid::~Grid()
{
    //dtor
}

Grid::Grid(const Grid& other)
{
    //copy ctor
}

Grid& Grid::operator=(const Grid& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void Grid::creategrid(const int &_nx, const int & _ny, const double & _dx, const double & _dy)
{
	geometric_parameters.nx = _nx;
	geometric_parameters.ny = _ny;
	geometric_parameters.dx = _dx;
	geometric_parameters.dy = _dy;
	cell_property.resize(geometric_parameters.nx);
	for (unsigned int ii = 0; ii < geometric_parameters.nx; ii++)
	{
		cell_property[ii].resize(geometric_parameters.ny);
	}
	for (unsigned int ii = 0; ii < geometric_parameters.nx; ii++)
	{
        for (unsigned int jj = 0; jj < geometric_parameters.ny; jj++)
        {
            cell_property[ii][jj].SetParent(this);
            cell_property[ii][jj].i = ii;
            cell_property[ii][jj].j = jj;
            cout<<ii<<", "<<jj<<", "<<cell_property[ii][jj].GetParent()<<endl;

        }
	}
}

void Grid::Show()
{
	for (unsigned int ii = 0; ii < geometric_parameters.nx; ii++)
	{
        for (unsigned int jj = 0; jj < geometric_parameters.ny; jj++)
        {
            cout<<cell_property[ii][jj].i<<", "<<cell_property[ii][jj].j<<", "<<cell_property[ii][jj].GetParent()<<endl;
        }
	}
}

void Grid::set_inverse_K_marginal_dist(int ninc)
{
	inverse_K_distribution.clear();
	double epsilon = 1e-8;

	for (int i = 0; i < ninc+1; i++)
	{
		double u = double(i) / double(ninc)*(1 - 2 * epsilon) + epsilon;
		if (i==0)
			inverse_K_distribution.append(u, inverse_cumulative_distribution(u));
		else
			inverse_K_distribution.append(u, inverse_cumulative_distribution(u));
	}

}

double Grid::inverse_cumulative_distribution(double u)
{
	double x2 = 1000;
	double x1 = 0.001;
	double tol = 1e-8;
	double err1 = cumulative_marginal_K_dist(x1) - u;
	double err2 = cumulative_marginal_K_dist(x2) - u;
	while (err1 > 0)
	{
		x1 /= 2;
		err1 = cumulative_marginal_K_dist(x1) - u;
	}

	while (err2 < 0)
	{
		x2 *= 2;
		err2 = cumulative_marginal_K_dist(x2) - u;
	}

	while (min(fabs(err1),fabs(err2)) > tol && fabs(x1-x2)>tol)
	{
		double slope = (err2 - err1) / (log(x2) - log(x1));
		double x_p = exp(log(x1) - err1 / slope);
		double err_p = cumulative_marginal_K_dist(x_p) - u;
		if (err_p > 0)
		{
			x2 = x_p;
			err2 = cumulative_marginal_K_dist(x2) - u;
		}
		else
		{
			x1 = x_p;
			err1 = cumulative_marginal_K_dist(x1) - u;
		}

	}
	if (fabs(err1) > fabs(err2))
		return x2;
	else
		return x1;
}

double Grid::cumulative_marginal_K_dist(double x)
{
	double out = 0;
	if (fieldgeneratorparameters.marginal_K_distribution_type == "lognormal")
	{
		int n = marginal_K_dist_params.size() / 3.0;

		for (int i = 0; i < n; i++)
			out += marginal_K_dist_params[3 * i] * 0.5*(1.0 + erf((log(x) - log(marginal_K_dist_params[3 * i + 1])) / (sqrt(2.0)*marginal_K_dist_params[3 * i + 2])));
	}
	return out;

}

void Grid::assign_K_gauss(int i, int j)
{
	CNormalDist ND;
	correl_mat_vec M = get_correll_matrix_vec(i, j);
	double mu;
	double sigma;
	if (M.V_RHS.num == 0)
	{
		mu = 0;
		sigma = 1;
	}
	else
	{
            CMatrix_arma M_inv = inv(M.M_22);
            mu = dotproduct(M_inv*M.V_21, M.V_RHS);
            sigma = 1.0 - dotproduct(M_inv*M.V_21, M.V_21);
	}

	double K_gauss = getnormalrand(mu, sigma);
	cell_property[i][j].k_det = true;
	n_k_dets++;
	cell_property[i][j].K_gauss = K_gauss;
	cell_property[i][j].K = map_to_KCDF(getnormalcdf(K_gauss));
}


void Grid::assign_K_gauss()
{
	set_inverse_K_marginal_dist(100);
	int n_filled = 0;
	srand(time(NULL));
	while (n_filled<geometric_parameters.nx*geometric_parameters.ny)
	{
        if (n_filled%100==0)
		set_progress_value(double(n_filled) / double(geometric_parameters.nx * geometric_parameters.ny));
            int i = unitrandom()*(geometric_parameters.nx-1) + 0.5;
            int j = unitrandom()*(geometric_parameters.ny-1) + 0.5;
            if (!cell_property[i][j].k_det)
            {
                assign_K_gauss(i, j);
                n_filled++;
            }
    }
}

correl_mat_vec Grid::get_correll_matrix_vec(int i, int j)
{
	correl_mat_vec M;
	vector<ijval> ij = get_closest_n_points(get_closest_K_dets(i, j, min(n_k_dets,fieldgeneratorparameters.max_number_of_points_used)+1), min(n_k_dets, fieldgeneratorparameters.max_number_of_points_used)+1);
#ifdef  USE_ARMA
	M.M_22 = CMatrix_arma(ij.size() - 1);
	M.V_21 = CVector_arma(ij.size() - 1);
	M.V_RHS = CVector_arma(ij.size() - 1);
#else
	M.M_22 = CMatrix(ij.size() - 1);
	M.V_21 = CVector(ij.size() - 1);
	M.V_RHS = CVector(ij.size() - 1);
#endif //  USE_ARMA
	for (int ii = 1; ii < int(ij.size()); ii++)
	{
		M.V_21[ii-1] = exp(-sqrt((i - ij[ii].i)*geometric_parameters.dx/ fieldgeneratorparameters.k_correlation_lenght_scale_x*(i - ij[ii].i)*geometric_parameters.dx/ fieldgeneratorparameters.k_correlation_lenght_scale_x + (j - ij[ii].j)*geometric_parameters.dy/ fieldgeneratorparameters.k_correlation_lenght_scale_y*(j - ij[ii].j)*geometric_parameters.dy/ fieldgeneratorparameters.k_correlation_lenght_scale_y) );
		M.V_RHS[ii - 1] = cell_property[ij[ii].i][ij[ii].j].K_gauss;
		for (int jj = 1; jj < int(ij.size()); jj++)
		{
#ifdef USE_ARMA
		M.M_22(ii-1,jj-1) = exp(-sqrt((ij[jj].i - ij[ii].i)*geometric_parameters.dx/ fieldgeneratorparameters.k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*geometric_parameters.dx/ fieldgeneratorparameters.k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*geometric_parameters.dy/ fieldgeneratorparameters.k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*geometric_parameters.dy/ fieldgeneratorparameters.k_correlation_lenght_scale_y));
#else
		M.M_22[ii - 1][jj - 1] = exp(-sqrt((ij[jj].i - ij[ii].i)*geometric_parameters.dx/ fieldgeneratorparameters.k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*geometric_parameters.dx/ fieldgeneratorparameters.k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*geometric_parameters.dy/ fieldgeneratorparameters.k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*geometric_parameters.dy/ fieldgeneratorparameters.k_correlation_lenght_scale_y));
#endif // USE_ARMA
		}
	}
	return M;
}

double Grid::map_to_KCDF(double u)
{
	return inverse_K_distribution.interpol(u);
}

void Grid::show_in_window(string s)
{
    #ifdef QT_version
    qDebug()<<QString::fromStdString(s);
    main_window->get_ui()->ShowOutput->append(QString::fromStdString(s));
    QApplication::processEvents();
    #else
    cout<<s<<endl;
    #endif // Qt_version
}

void Grid::set_progress_value(double s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s*100 << "%                                     ";
}

void Grid::set_progress_value(string s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r" << s << "                                     ";
}

vector<ijval> Grid::get_closest_n_points(vector<ijval> vec, int n)
{
	vector<ijval> out;
	vector<bool> extracted(vec.size());
	for (int i = 0; i < int(vec.size()); i++) extracted[i] = false;
	int smallest_dist = -1;

	for (int i = 0; i < n; i++)
	{
		double min_dist = 1e12;
		for (int j = 0; j < int(vec.size()); j++)
			if ((vec[j].val < min_dist) && (extracted[j]==false))
			{
				smallest_dist = j;
				min_dist = vec[j].val;
			}
		out.push_back(vec[smallest_dist]);
		extracted[smallest_dist] = true;
	}

	return out;
}

vector<ijval> Grid::get_closest_K_dets(int i, int j, int n)
{
	vector<ijval> out;
	double max_dist = 0;
	for (int k = 1; k < max(geometric_parameters.nx, geometric_parameters.ny); k++)
	{
		double min_dist = 1e12;
		int jj = j - k;
		if (jj>=0)
			for (int ii = max(i - k,0); ii <= min(i + k,geometric_parameters.nx-1); ii++)
			{
				double dist2 = ((i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x*(i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x + (j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y*(j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (cell_property[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		jj = j + k;
		if (jj < geometric_parameters.ny)
			for (int ii = max(i - k, 0); ii <= min(i + k, geometric_parameters.nx - 1); ii++)
			{
				double dist2 = ((i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x*(i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x + (j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y*(j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (cell_property[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		int ii = i - k;
		if (ii >= 0)
			for (int jj = max(j - k+1, 0); jj <= min(j + k-1, geometric_parameters.ny - 1); jj++)
			{
				double dist2 = ((i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x*(i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x + (j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y*(j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (cell_property[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		ii = i + k;
		if (ii < geometric_parameters.nx)
			for (int jj = max(j - k + 1, 0); jj <= min(j + k - 1, geometric_parameters.ny - 1); jj++)
			{
				double dist2 = ((i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x*(i - ii)*geometric_parameters.dx/fieldgeneratorparameters.k_correlation_lenght_scale_x + (j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y*(j - jj)*geometric_parameters.dy/fieldgeneratorparameters.k_correlation_lenght_scale_y);
				min_dist = min(min_dist, dist2);
				if (cell_property[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
                                            ijval pp;
                                            pp.i = ii;
                                            pp.j = jj;
                                            pp.val = sqrt(dist2);
                                            out.push_back(pp);
                                            max_dist = max(dist2, max_dist);
					}
				}

			}
		if ((int(out.size()) >= n) && min_dist > max_dist)
		{
                    k = max(geometric_parameters.nx, geometric_parameters.ny);
                    ijval pp;
                    pp.i = i;
                    pp.j = j;
                    pp.val = 0;
                    out.push_back(pp);
                    return out;
		}
	}
	ijval pp;
	pp.i = i;
	pp.j = j;
	pp.val = 0;
	out.push_back(pp);
	return out;
}

double Grid::max_K()
{
	double max_k = -1e23;
	for (int i = 0; i < geometric_parameters.nx; i++)
		for (int j = 0; j < geometric_parameters.ny; j++)
			max_k = max(max_k, cell_property[i][j].K);

	return max_k;

}
double Grid::min_K()
{
	double min_k = 1e23;
	for (int i = 0; i < geometric_parameters.nx; i++)
		for (int j = 0; j < geometric_parameters.ny; j++)
			min_k = min(min_k, cell_property[i][j].K);
	return min_k;
}

void Grid::renormalize_k()
{
	CBTC K = get_marginal_K_in_array();
	double mu = K.mean();
	double std = K.std();
	for (int i = 0; i < geometric_parameters.nx; i++)
		for (int j = 0; j < geometric_parameters.ny; j++)
			cell_property[i][j].K_gauss = (cell_property[i][j].K_gauss - mu) / std;
	remap_K();
}

void Grid::remap_K()
{
	for (int i = 0; i < geometric_parameters.nx; i++)
		for (int j = 0; j < geometric_parameters.ny; j++)
			cell_property[i][j].K = map_to_KCDF(getnormalcdf(cell_property[i][j].K_gauss));

}

CTimeSeries Grid::get_marginal_K_in_array()
{
	CBTC out;
	for (int i = 0; i < geometric_parameters.nx; i++)
		for (int j = 0; j < geometric_parameters.ny; j++)
			out.append(i + j * 10000, cell_property[i][j].K_gauss);

	return out;
}

void Grid::AssignValues(const string &option, vector<double> &parameters)
{
    if (option=="linear_in_x")
    {
        for (int i = 0; i < geometric_parameters.nx; i++)
            for (int j = 0; j < geometric_parameters.ny; j++)
                cell_property[i][j].value = parameters[0]*double(i)/double(geometric_parameters.nx);
    }

}

Cell_Property* Grid::Cell(Position *p)
{
     return &cell_property[p->I()][p->J()];
}

CTimeSeries Grid::GetValueAsArray()
{
    CTimeSeries valuearray;
    for (int i = 0; i < geometric_parameters.nx; i++)
        for (int j = 0; j < geometric_parameters.ny; j++)
            valuearray.append(j+10000*i,cell_property[i][j].value);
    return valuearray;

}
CTimeSeries Grid::GetValueDistribution(unsigned int number_of_bins)
{
    return GetValueAsArray().distribution(number_of_bins);
}
void Grid::SetValueRanks(unsigned int number_of_bins)
{
    CTimeSeries inverse_cummulative = GetValueDistribution(number_of_bins).getcummulative().make_uniform(number_of_bins);
    for (int i = 0; i < geometric_parameters.nx; i++)
        for (int j = 0; j < geometric_parameters.ny; j++)
            cell_property[i][j].ValueRank = inverse_cummulative.interpol(cell_property[i][j].value);
}

CMatrix Grid::GetMatrix(const string &property)
{
    CMatrix valuearray(geometric_parameters.nx,geometric_parameters.ny);
    for (unsigned int i = 0; i < geometric_parameters.nx; i++)
        for (unsigned int j = 0; j < geometric_parameters.ny; j++)
            valuearray(i,j) = cell_property[i][j].PropertyValue(property);

    return valuearray;
}



