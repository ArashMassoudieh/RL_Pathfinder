#include "Pathway.h"
#include "Grid.h"
#include "Utilities.h"
#include "ReinforcementLearner.h"
#include <iostream>
#include "PathwaySet.h"

using namespace std;

Pathway::Pathway(Grid *_parent)
{
    parent = _parent;
}

Pathway::~Pathway()
{
    //dtor
}

Pathway::Pathway(const Pathway& other)
{
    points = other.points;
}

Pathway& Pathway::operator=(const Pathway& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    points = rhs.points;
    return *this;
}

Position& Pathway::OneStepMove(bool store_results, bool writeout)
{
    if (!Done)
    {
        if (CurrentPoint().I()==grid_nx()-1)
        {
            Done = true;
            if (writeout) cout<<"Done!"<<endl;
        }
        CVector values = ValueofActions();
        parent->Cell(&Current_Point)->SetValue(GetReinforcementLearner()->Parameters.Gamma*values.max()+parent->Cell(&Current_Point)->Reward());
        if (unitrandom()>GetReinforcementLearner()->Parameters.epsilon)
        {
            int actionselected = aquiutils::randompick(values.maxelements());
            Current_Point = ActionToPosition(actionselected);
            if (store_results)
                points.push_back(Current_Point);
        }
        else
        {
            int actionselected = rand()%values.num;
            Current_Point = ActionToPosition(actionselected);
            if (store_results)
                points.push_back(Current_Point);
        }

        return Current_Point;
    }
}

void Pathway::Initialize(const string s, unsigned int initial_i)
{
    //srand (time(NULL));
    if (parent==nullptr)
    {
        cout<<"The pathways lacks a parent!"<<endl;
        return;
    }
    Done = false;
    Current_Point.SetPathway(this);
    double u1 = RandomNumberGenerator()->unitrandom();
    double u2 = RandomNumberGenerator()->unitrandom();
    if (s=="random_x")
        Current_Point.SetI(int(u1*parent->Geometric_Parameters().nx));
    else
        Current_Point.SetI(initial_i);
    Current_Point.SetJ(int(u2*parent->Geometric_Parameters().ny));
    Current_Point.SetX(Current_Point.I()*grid_dx()+grid_dx()*RandomNumberGenerator()->unitrandom()*0.05);
    Current_Point.SetY(Current_Point.J()*grid_dy()+grid_dy()*RandomNumberGenerator()->unitrandom()*0.05);
    points.push_back(Current_Point);
}

CVector Pathway::ValueofActions()
{
    CVector valuesofactions(pathwayset->numberofdirections);
    for (unsigned int i=0; i<pathwayset->numberofdirections; i++)
    {
        valuesofactions[i] = value(ActionToPosition(i));
    }
    return valuesofactions;
}

Position Pathway::ActionToPosition(unsigned int a)
{
    Position p = Current_Point;
    if (a==0)
        p.SetI(min(grid_nx()-1,p.I()+1));
    if (a==1)
        p.SetI(max(0,p.I()-1));
    if (a==2)
        p.SetJ(min(grid_ny()-1, p.J()+1));
    if (a==3)
        p.SetJ(max(0, p.J()-1));
    if (a==4)
    {
        double a = 1;//do nothing
    }

    return p;
}

double& Pathway::value(int i, int j)
{
    return parent->Cell(i,j).value;

}

double& Pathway::value(const Position &p)
{
    return parent->Cell(p.I(),p.J()).value;
}

double Pathway::grid_dx() {
    if (parent != nullptr)
        return parent->Geometric_Parameters().dx;
    else
        return 0;
}
double Pathway::grid_dy() {
    if (parent != nullptr)
        return parent->Geometric_Parameters().dy;
    else
        return 0;
}
int Pathway::grid_nx() {
    if (parent != nullptr)
        return parent->Geometric_Parameters().nx;
    else
        return 0;
}
int Pathway::grid_ny() {
    if (parent != nullptr)
        return parent->Geometric_Parameters().ny;
    else
        return 0;
}

ReinforcementLearner *Pathway::GetReinforcementLearner()
{
    if (parent != nullptr)
        return parent->GetReinforcementLearner();
    else
        return nullptr;
}

Grid* Pathway::grid()
{
    return parent;
}

void Pathway::SetGrid(Grid *g)
{
    parent = g;
}

void Pathway::SetPathwaySet(PathwaySet *g)
{
    pathwayset = g;
}

void Pathway::Clear()
{
    points.clear();
}

CNormalDist* Pathway::RandomNumberGenerator() {
    return &parent->RandomNumberGenerator;
}
