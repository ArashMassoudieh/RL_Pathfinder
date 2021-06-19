#include "Cell_Property.h"
#include <iostream>
#include "Grid.h"

using namespace std;

Cell_Property::Cell_Property(Grid *_parent)
{
    parent = _parent;
}

Cell_Property::~Cell_Property()
{
    //dtor
}

Cell_Property::Cell_Property(const Cell_Property& other)
{
    K = other.K;
    K_gauss = other.K_gauss;
    k_det = other.k_det;
    value = other.value;
    ValueRank = other.ValueRank;
    i = other.i;
    j = other.j;
}

Cell_Property& Cell_Property::operator=(const Cell_Property& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    K = rhs.K;
    K_gauss = rhs.K_gauss;
    k_det = rhs.k_det;
    value = rhs.value;
    i = rhs.i;
    j = rhs.j;
    ValueRank = rhs.ValueRank;
    return *this;
}

double Cell_Property::PropertyValue(const string quantity)
{
    if (quantity == "K")
        return K;
    if (quantity == "K_gauss")
        return K_gauss;
    if (quantity == "value")
        return value;
    if (quantity == "reward")
        return Reward();
    if (quantity == "valuerank")
        return ValueRank;

}

void Cell_Property::SetValue(double val)
{
    value = val;
}

double Cell_Property::Reward()
{
    double reward = -K;
    //cout<<"Reward:" <<i<<","<<j<<","<<parent<<endl;
    if (i==parent->Geometric_Parameters().nx-1)
        reward += parent->GetReinforcementLearner()->Parameters.UltimateReward;
    return reward;
}

void Cell_Property::SetParent(Grid *_grid)
{
    parent = _grid;
}
