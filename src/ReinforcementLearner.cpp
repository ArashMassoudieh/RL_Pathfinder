#include "ReinforcementLearner.h"

ReinforcementLearner::ReinforcementLearner()
{
    //ctor
}

ReinforcementLearner::~ReinforcementLearner()
{
    //dtor
}

ReinforcementLearner::ReinforcementLearner(const ReinforcementLearner& other)
{
    Parameters = other.Parameters;
}

ReinforcementLearner& ReinforcementLearner::operator=(const ReinforcementLearner& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    Parameters = rhs.Parameters;
    return *this;
}
