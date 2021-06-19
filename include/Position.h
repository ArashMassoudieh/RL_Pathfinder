#ifndef POSITION_H
#define POSITION_H
#include <vector>
using namespace std;

class Pathway;
class Grid;

class Position
{
    public:
        Position(Pathway* pthway=nullptr);
        virtual ~Position();
        Position(const Position& other);
        Position& operator=(const Position& other);
        double X() const;
        double Y() const;
        int I() const {return i;}
        int J() const {return j;}
        void SetX(double xx) {x = xx;}
        void SetY(double yy) {y = yy;}
        void SetI(int ii) {i = ii;}
        void SetJ(int jj) {j = jj;}
        Grid* grid();
        Pathway* pathway();
        void SetPathway(Pathway *pthway);
    protected:

    private:
        double x;
        double y;
        int i;
        int j;
        vector<double> values;
        Pathway *parent = nullptr;

};

#endif // POSITION_H
