#include <iostream>
#include "Grid.h"
#include "Pathway.h"
#include "PathwaySet.h"


using namespace std;

int main()
{
    int number_of_paths = 100000;
    int number_of_epochs = 200;
    int number_of_steps = 10000;
    int nx = 50;
    int ny = 50;
    Grid G;
    G.creategrid(nx,ny,1/double(nx),1/double(ny));
    G.fieldgeneratorparameters.k_correlation_lenght_scale_x=0.5;
    G.fieldgeneratorparameters.k_correlation_lenght_scale_y=0.5;
    G.fieldgeneratorparameters.max_number_of_points_used = 10;
    G.fieldgeneratorparameters.marginal_K_distribution_type = "lognormal";
    vector<double> distribution_parameters;
    distribution_parameters.push_back(0.9);
    distribution_parameters.push_back(10);
    distribution_parameters.push_back(1);
    distribution_parameters.push_back(0.1);
    distribution_parameters.push_back(0.1);
    distribution_parameters.push_back(1);
    G.Set_Marginal_K_Distribution_Parameters(distribution_parameters);
    G.set_inverse_K_marginal_dist(100);
    G.assign_K_gauss();
    G.renormalize_k();
    G.write_K_field_to_vtp("/mnt/3rd900/Projects/RL_path_results/results/test.vtp",0.0);
    vector<double> value_parameters;
    value_parameters.push_back(0.1);
    ReinforcementLearner RL;
    RL.Parameters.epsilon = 0.1;
    RL.Parameters.UltimateReward = 1*exp(pow(1,2)/2)*double(nx);
    G.SetReinforcementLearner(&RL);
    //G.AssignValues("linear_in_x",value_parameters);

    G.write_field_to_vtp("reward","/mnt/3rd900/Projects/RL_path_results/results/rewards.vtp",0);
    PathwaySet pthwyset(&G);
    pthwyset.numberofdirections = 5;
    for (unsigned int j=0; j<number_of_epochs; j++)
    {
        int numberofsteps=0;
        pthwyset.Initialize(number_of_paths, "random_x");
        unsigned int numberofdones = 0;
        for (unsigned int i=0; i<number_of_steps; i++)
        {
            numberofdones = pthwyset.OneStepMove(false, false);
            numberofsteps += number_of_paths-numberofdones;
            if (numberofdones == number_of_paths)
                break;
        }
        cout<<j<<":"<<numberofsteps<<","<<numberofdones<<endl;;
        G.GetMatrix("value").writetofile("/mnt/3rd900/Projects/RL_path_results/results/value_matrix_" + aquiutils::numbertostring(j)+".txt");
        G.GetValueDistribution(100).writefile("/mnt/3rd900/Projects/RL_path_results/results/Value_Distribution_"+aquiutils::numbertostring(j)+".txt");
        G.SetValueRanks(50);
        G.GetMatrix("valuerank").writetofile("/mnt/3rd900/Projects/RL_path_results/results/rank_matrix_" + aquiutils::numbertostring(j)+".txt");
        G.write_field_to_vtp("value","/mnt/3rd900/Projects/RL_path_results/results/values_"+aquiutils::numbertostring(j)+".vtp",0);
        G.write_field_to_vtp("valuerank","/mnt/3rd900/Projects/RL_path_results/results/ranks_"+aquiutils::numbertostring(j)+".vtp",0);
    }

    //pthwyset.SaveToVTP("path.vtp",1,0,1);

    G.write_field_to_vtp("value","/mnt/3rd900/Projects/RL_path_results/results/values.vtp",0);

    PathwaySet pathset(&G);
    pathset.Initialize(100,"none_random_x", 0);
    pathset.numberofdirections = 4;
    RL.Parameters.epsilon = 0;
    unsigned int numberofdones = 0;
    while(numberofdones<100)
    {
        numberofdones = pathset.OneStepMove(true,false);
        cout<<"n_done="<<numberofdones<<endl;
    }

    pathset.SaveToVTP("/mnt/3rd900/Projects/RL_path_results/results/final_path.vtp",1,0,1);

}
