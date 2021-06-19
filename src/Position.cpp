#include "Position.h"
#include "Pathway.h"
#include "Grid.h"

Position::Position(Pathway* pthway)
{
    parent = pthway;
}

Position::~Position()
{
    //dtor
}

Position::Position(const Position& other)
{
    x = other.x;
    y = other.y;
    i = other.i;
    j = other.j;
    parent = other.parent;
    values = other.values;

}

Position& Position::operator=(const Position& other)
{
    if (this == &other) return *this; // handle self assignment
    x = other.x;
    y = other.y;
    i = other.i;
    j = other.j;
    parent = other.parent;
    values = other.values;
    return *this;
}

Grid* Position::grid()
{
    return parent->grid();
}

Pathway* Position::pathway()
{
    return parent;
}

void Position::SetPathway(Pathway *pthway)
{
    parent = pthway;
}

double Position::X() const
{
    if (parent!=nullptr)
        return parent->grid_dx()*i;
}

double Position::Y() const
{
    if (parent!=nullptr)
        return parent->grid_dy()*j;
}
