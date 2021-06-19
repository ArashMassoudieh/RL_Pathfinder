#ifndef GRID_H
#define GRID_H

#include <vector>
#include "Matrix_arma.h"
#include "Vector_arma.h"
#include "BTCSet.h"
#include "vtk.h"
#include "ReinforcementLearner.h"
#include "Cell_Property.h"
#include "Position.h"

struct ijval
{
	int i;
	int j;
	double val;
};

struct correl_mat_vec
{
#ifdef  USE_ARMA

	CMatrix_arma M_22;
	CVector_arma V_21;
	CVector_arma V_RHS;
#else
	CMatrix M_22;
	CVector V_21;
	CVector V_RHS;
#endif //  arma
};

struct Grid_Geometric_Parameters
{
	int nx, ny;
	double dx, dy;
	int nx_data, ny_data;
};

struct Field_Generator_Parameters
{
	int max_number_of_points_used;
	double k_correlation_lenght_scale_x;
	double k_correlation_lenght_scale_y;
	string marginal_K_distribution_type;
};

using namespace std;
class Grid
{
    public:
        Grid();
        virtual ~Grid();
        Grid(const Grid& other);
        Grid& operator=(const Grid& other);
        void creategrid(const int &_nx, const int & _ny, const double & _dx, const double & _dy);
        void Show();
        void set_inverse_K_marginal_dist(int ninc);
        Field_Generator_Parameters fieldgeneratorparameters;
        void assign_K_gauss();
        void show_in_window(string s);
        void set_progress_value(double s);
        void set_progress_value(string s);
        void Set_Marginal_K_Distribution_Parameters(const vector<double> &parameters) {marginal_K_dist_params = parameters;}
        void write_K_field_to_vtp(string filename, double z_scale=1);
        void write_field_to_vtp(const string &quantity, const string &filename, double z_factor);
        void renormalize_k();
        void remap_K();
        Grid_Geometric_Parameters &Geometric_Parameters() {return geometric_parameters;}
        Cell_Property &Cell(unsigned int i, unsigned int j) {return cell_property[i][j];}
        Cell_Property* Cell(Position *p);
        void AssignValues(const string &option, vector<double> &parameters);
        void SetReinforcementLearner(ReinforcementLearner *rl) {reinforcementlearner = rl;}
        ReinforcementLearner* GetReinforcementLearner() {return reinforcementlearner;}
        CMatrix GetMatrix(const string &property);
        CTimeSeries GetValueAsArray();
        CTimeSeries GetValueDistribution(unsigned int number_of_bins);
        void SetValueRanks(unsigned int number_of_bins);
        CNormalDist RandomNumberGenerator;
    protected:

    private:
        Grid_Geometric_Parameters geometric_parameters;
        vector<ijval> get_closest_K_dets(int i, int j, int n);
        correl_mat_vec get_correll_matrix_vec(int i, int j);
        void assign_K_gauss(int i, int j);
        vector<vector<Cell_Property>> cell_property;
        vector<double> marginal_K_dist_params;
        double inverse_cumulative_distribution(double u);
        double cumulative_marginal_K_dist(double x);
        CTimeSeries inverse_K_distribution;
        double map_to_KCDF(double u);
        int n_k_dets = 0;
        vector<ijval> get_closest_n_points(vector<ijval> vec, int n);
        CTimeSeries get_marginal_K_in_array();
        double max_K();
        double min_K();
        ReinforcementLearner *reinforcementlearner;

};

#endif // GRID_H
