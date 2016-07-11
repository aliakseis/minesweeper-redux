// MineSweeperRedux.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>

#include <set>
#include <vector>

#include <algorithm>
#include <numeric>

using std::cout;

using std::set;
using std::vector;
using std::min;
using std::max;

using std::sort;

const char squares[10][11] =
{
    "*****2*4*3",
    "3**0******",
    "**101***3*",
    "3****1***2",
    "2******2**",
    "1*5*2*0**1",
    "******1***",
    "****1**3**",
    "*23**1***2",
    "1*3*****2*",
};

struct Cell
{
    const int x, y;
    Cell(int x_, int y_) : x(x_), y(y_) {}
    bool operator < (const Cell& other) const
    {
        return x < other.x || x == other.x && y < other.y;
    }
};

struct Group
{
    set<Cell> neigbors;
    int countMin, countMax;

    bool test;

    Group() : test(false) {}

    bool operator < (const Group& other) const
    {
        return neigbors.size() < other.neigbors.size();
    }
};

bool HandleUnambiguous(vector<Group>& groups, set<Cell>& result)
{
    vector<Group>::iterator itEnd = groups.end();
    for (vector<Group>::iterator it = groups.begin()
        ; it != itEnd
        ; ++it)
    {
        const bool areZeros = 0 == it->countMax;
        if (areZeros || it->neigbors.size() == it->countMin)
        {
            set<Cell> cells;
            cells.swap(it->neigbors);
            groups.erase(it);

            if (!areZeros)
                result.insert(cells.begin(), cells.end());

            for (int i = groups.size(); --i >= 0;)
            {
                set<Cell>::iterator itEnd = cells.end();
                for (set<Cell>::iterator it = cells.begin()
                        ; it != itEnd
                        ; ++it)
                    if (groups[i].neigbors.erase(*it) && !areZeros)
                    {
                        groups[i].countMin--;
                        groups[i].countMax--;
                    }

                if (0 == groups[i].neigbors.size())
                    groups.erase(groups.begin() + i);
            }

            return true;
        }
    }
    return false;
}

/*
bool HandleCommon(vector<Group>& groups)
{
    sort(groups.begin(), groups.end());
    for (int i = groups.size(); --i > 0; )
        for (int j = i; --j >= 0; )
        {
            set<Cell> test = groups[i].neigbors;
            test.insert(groups[j].neigbors.begin(), groups[j].neigbors.end());
            if (test.size() != groups[i].neigbors.size())
                continue;

            set<Cell>::iterator itEnd = groups[j].neigbors.end();
            for (set<Cell>::iterator it = groups[j].neigbors.begin()
                    ; it != itEnd
                    ; ++it)
                groups[i].neigbors.erase(*it);

            groups[i].count -= groups[j].count;

            groups.erase(groups.begin() + j);

            return true;
        }

    return false;
}
*/

bool HandleIntersections(vector<Group>& groups)
{
    bool ok = false;

    for (int i = groups.size(); --i > 0; )
        for (int j = i; --j >= 0; )
        {
            if (groups[i].test && groups[j].test)
                __asm int 3;

            Group* pI = &groups[i];
            Group* pJ = &groups[j];

            set<Cell> setI, setJ, intersect;
            set<Cell>::iterator itI(groups[i].neigbors.begin()), itJ(groups[j].neigbors.begin());
            while (itI != groups[i].neigbors.end() && itJ != groups[j].neigbors.end())
                if (*itI < *itJ)
                {
                    setI.insert(*itI);
                    ++itI;
                }
                else if (*itJ < *itI)
                {
                    setJ.insert(*itJ);
                    ++itJ;
                }
                else
                {
                    intersect.insert(*itI);
                    ++itI;
                    ++itJ;
                }

            if (intersect.size() < 2)
                continue;

            ok = true;

            setI.insert(itI, groups[i].neigbors.end());
            setJ.insert(itJ, groups[j].neigbors.end());

            int intersectMax = min((int)intersect.size(), min(groups[i].countMax, groups[j].countMax));

            int intersectMin = max(0, 
                (int)(intersect.size() - min(groups[i].neigbors.size() - groups[i].countMin, groups[j].neigbors.size() - groups[j].countMin)));

            {
                Group groupI;
                groupI.neigbors.swap(setI);
                groupI.countMin = groups[i].countMin - intersectMax;
                groupI.countMax = groups[i].countMax - intersectMin;
                groups.push_back(groupI);
            }
            {
                Group groupJ;
                groupJ.neigbors.swap(setJ);
                groupJ.countMin = groups[j].countMin - intersectMax;
                groupJ.countMax = groups[j].countMax - intersectMin;
                groups.push_back(groupJ);
            }
            {
                Group group;
                group.neigbors.swap(intersect);
                group.countMin = intersectMin;
                group.countMax = intersectMax;
                groups.push_back(group);
            }
        }

    return ok;
}

int _tmain(int argc, _TCHAR* argv[])
{
    vector<Group> groups;
    set<Cell> result;

    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
        {
            char c = squares[i][j];
            if (c >= '0' && c <= '8')
            {
                Group group;
                group.test = (i == 3 && j == 0 || i == 2 && j == 2);
                group.countMin = group.countMax = c - '0';
                for (int i1 = max(i-1, 0); i1 <= min(i+1, 9); ++i1)
                    for (int j1 = max(j-1, 0); j1 <= min(j+1, 9); ++j1)
                        if ((i != i1 || j != j1) && (squares[i1][j1] < '0' || squares[i1][j1] > '8'))
                            group.neigbors.insert(Cell(i1, j1));

                groups.push_back(group);
            }
        }

    while (HandleUnambiguous(groups, result))
        ;

    HandleIntersections(groups);

    while (HandleUnambiguous(groups, result))
        ;

    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            cout << ((result.find(Cell(i, j)) != result.end()) ? 'M' : squares[i][j]);
        }
        cout << '\n';
    }

	return 0;
}

