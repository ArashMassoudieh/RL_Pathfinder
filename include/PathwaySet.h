#ifndef PATHWAYSET_H
#define PATHWAYSET_H

#include <vector>
#include "Pathway.h"

class Grid;

using namespace std;
class PathwaySet
{
    public:
        PathwaySet(Grid *_parent);
        virtual ~PathwaySet();
        PathwaySet(const PathwaySet& other);
        PathwaySet& operator=(const PathwaySet& other);
        Pathway *operator[](int i);
        void SaveToVTP(vtkSmartPointer<vtkPolyDataMapper> mapper, string filename);
        vtkSmartPointer<vtkPolyDataMapper> ConvertToVTP(double z_factor, double offset);
        void SaveToVTP(string filename, double z_factor, double offset, int interval);
        void Append(Pathway &p);
        void Initialize(int number_of_paths, const string &s = "random_x", int starting_point = 0);
        unsigned int OneStepMove(bool store_results = true, bool periodic=true);
        void Clear() {paths.clear();};
        unsigned int numberofdirections = 5;
        unsigned int number_reached_goal=0;
    protected:

    private:
        vector<Pathway> paths;
        Grid *parent;

};

#endif // PATHWAYSET_H
