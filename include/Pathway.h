#ifndef PATHWAY_H
#define PATHWAY_H
#include "Position.h"
#include "NormalDist.h"
#include "vtk.h"

class Grid;
class ReinforcementLearner;
class PathwaySet;

class Pathway
{
    public:
        Pathway(Grid* _parent = nullptr);
        virtual ~Pathway();
        Pathway(const Pathway& other);
        Pathway& operator=(const Pathway& other);
        Position& OneStepMove(bool store_results=false, bool writedone = false);
        void Initialize(const string s = "random_x", unsigned int i=0);
        double grid_dx();
        double grid_dy();
        int grid_nx();
        int grid_ny();
        int PickAction();
        CVector ValueofActions();
        Position ActionToPosition(unsigned int a);
        double& value(int i, int j);
        double& value(const Position &p);
        ReinforcementLearner *GetReinforcementLearner();
        vtkSmartPointer<vtkPolyData> ToVTP(double z_factor, double offset, int ID);
        Grid* grid();
        void SetGrid(Grid *g);
        void SetPathwaySet(PathwaySet *g);
        Position &CurrentPoint() {return Current_Point;}
        bool Done = false;
        void Clear();
    protected:

    private:
        vector<Position> points;
        Grid* parent;
        Position Current_Point;
        CNormalDist* RandomNumberGenerator();
        PathwaySet *pathwayset;
};

#endif // PATHWAY_H
