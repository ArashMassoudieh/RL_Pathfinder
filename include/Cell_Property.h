#ifndef CELL_PROPERTY_H
#define CELL_PROPERTY_H
#include <string>

class Grid;

using namespace std;
class Cell_Property
{
    public:
        Cell_Property(Grid *_parent = nullptr);
        virtual ~Cell_Property();
        Cell_Property(const Cell_Property& other);
        Cell_Property& operator=(const Cell_Property& other);
        double K;
        double K_gauss;
        bool k_det = false;
        double value = 0;
        double PropertyValue(const string quantity);
        void SetValue(double val);
        double ValueRank;
        double Reward();
        unsigned int i;
        unsigned int j;
        void SetParent(Grid *_grid);
        Grid* GetParent() {return parent;}

    protected:

    private:
        Grid *parent = nullptr;
};

#endif // CELL_PROPERTY_H
