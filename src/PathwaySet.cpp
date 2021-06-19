#include "PathwaySet.h"
#include "Grid.h"

PathwaySet::PathwaySet(Grid* _parent)
{
    parent = _parent;
}

PathwaySet::~PathwaySet()
{
    //dtor
}

PathwaySet::PathwaySet(const PathwaySet& other)
{
    paths = other.paths;
}

PathwaySet& PathwaySet::operator=(const PathwaySet& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    paths = rhs.paths;
    return *this;
}

void PathwaySet::Append(Pathway &p)
{
    paths.push_back(p);
    paths[paths.size()-1].SetGrid(parent);
    paths[paths.size()-1].SetPathwaySet(this);
}

void PathwaySet::Initialize(int number_of_paths, const string &s, int starting_point)
{
    number_reached_goal = 0;
    paths.resize(number_of_paths);
    for (unsigned int i=0; i<number_of_paths; i++)
    {
        paths[i].SetGrid(parent);
        paths[i].Clear();
        paths[i].SetPathwaySet(this);
        paths[i].Initialize(s, starting_point);
    }
}

unsigned int PathwaySet::OneStepMove(bool store_results, bool periodic)
{
    if (periodic)
    {
        for (unsigned int i=0; i<paths.size(); i++)
        {
            if (!paths[i].Done)
                paths[i].OneStepMove(store_results);
            else
            {
                number_reached_goal++;
                paths[i].Clear();
                paths[i].Initialize(0);
                paths[i].Done=false;
            }
        }
    }
    else
    {
        number_reached_goal = 0;
        for (unsigned int i=0; i<paths.size(); i++)
        {
            if (!paths[i].Done)
                paths[i].OneStepMove(store_results);
            else
            {
                number_reached_goal++;
            }
        }
    }
    return number_reached_goal;
}


