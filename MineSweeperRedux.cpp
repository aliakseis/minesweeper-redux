// MineSweeperRedux.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>

#include <set>
#include <vector>

#include <algorithm>
#include <numeric>

#include <shlwapi.h>
#pragma comment(lib, "shlwapi")

#include <assert.h>

using std::cout;
using std::endl;

using std::pair;
using std::set;
using std::vector;
using std::min;
using std::max;

using std::sort;

using std::swap;

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
    //const 
        int x, y;
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

    bool operator < (const Group& other) const
    {
        return neigbors.size() < other.neigbors.size()
            || neigbors.size() == other.neigbors.size()
            && lexicographical_compare(neigbors.begin(), neigbors.end(), other.neigbors.begin(), other.neigbors.end());
    }
};


class Prestidigitator
{
public:
    vector<set<Cell> > results;

    void ProduceResults(vector<Group>& groups)
    {
        set<Cell> result;
        ProduceResult(groups, result);
    }

private:
    template<typename T>
    bool Reduce(vector<Group>& groups, T begin, T end, const bool areZeros)
    {
        for (int i = (int)groups.size(); --i >= 0;)
        {
            for (T it = begin; it != end; ++it)
                if (groups[i].neigbors.erase(*it))
                    if (!areZeros)
                    {
                        if (groups[i].countMin > 0)
                            groups[i].countMin--;
                        groups[i].countMax--;
                    }
                    else
                    {
                        if (groups[i].countMax > (int)groups[i].neigbors.size())
                            groups[i].countMax = (int)groups[i].neigbors.size();
                    }

            if (groups[i].countMin > groups[i].countMax)
                return false;

            if (0 == groups[i].neigbors.size())
                groups.erase(groups.begin() + i);
        }

        return true;
    }

    bool Reduce(vector<Group>& groups, const set<Cell>& cells, const bool areZeros)
    {
        return Reduce(groups, cells.begin(), cells.end(), areZeros);
    }

    bool HandleUnambiguous(vector<Group>& groups, set<Cell>& result)
    {
        bool stop;
        do
        {
            stop = true;
            vector<Group>::iterator itEnd = groups.end();
            for (vector<Group>::iterator it = groups.begin()
                ; it != itEnd
                ; ++it)
            {
                if (it->countMax < 0 || it->countMin > (int)it->neigbors.size())
                    return false;
                const bool areZeros = 0 == it->countMax;
                if (areZeros || it->neigbors.size() == it->countMin)
                {
                    set<Cell> cells;
                    cells.swap(it->neigbors);
                    groups.erase(it);

                    if (!areZeros)
                        result.insert(cells.begin(), cells.end());

                    if (!Reduce(groups, cells, areZeros))
                        return false;

                    stop = false;
                    break;
                }
            }
        }
        while (!stop);

        return true;
    }

    bool Update(set<Group>& groups, const Group& newGroup)
    {
        pair<set<Group>::iterator, bool> inserted = groups.insert(newGroup);
        if (inserted.second)
            return true;

        bool result = false;
        if (inserted.first->countMin < newGroup.countMin)
        {
            result = true;
            inserted.first->countMin = newGroup.countMin;
        }
        if (inserted.first->countMax > newGroup.countMax)
        {
            result = true;
            inserted.first->countMax = newGroup.countMax;
        }

        return result;
    }


    void DoHandleIntersections(set<Group>::iterator it1, set<Group>::iterator it2, set<Group>& newGroups)
    {
        set<Cell> setI, setJ, intersect;
        set<Cell>::iterator itI(it1->neigbors.begin()), itJ(it2->neigbors.begin());
        while (itI != it1->neigbors.end() && itJ != it2->neigbors.end())
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

        if (intersect.size() < 2 
                || it1->neigbors.size() - intersect.size() < 2
                && it2->neigbors.size() - intersect.size() < 2)
            return;

        setI.insert(itI, it1->neigbors.end());
        setJ.insert(itJ, it2->neigbors.end());

        int intersectMax = min((int)intersect.size(), min(it1->countMax, it2->countMax));

        int intersectMin = max(0, 
            (int)(intersect.size() - min(it1->neigbors.size() - it1->countMin, it2->neigbors.size() - it2->countMin)));

        if (setI.size() > 1)
        {
            Group groupI;
            groupI.neigbors.swap(setI);
            groupI.countMin = max(it1->countMin - intersectMax, 0);
            groupI.countMax = min(it1->countMax - intersectMin, (int)groupI.neigbors.size());
            Update(newGroups, groupI);
        }
        if (setJ.size() > 1)
        {
            Group groupJ;
            groupJ.neigbors.swap(setJ);
            groupJ.countMin = max(it2->countMin - intersectMax, 0);
            groupJ.countMax = min(it2->countMax - intersectMin, (int)groupJ.neigbors.size());
            Update(newGroups, groupJ);
        }
        {
            Group group;
            group.neigbors.swap(intersect);
            group.countMin = intersectMin;
            group.countMax = intersectMax;
            Update(newGroups, group);
        }
    }

    void DoHandleIntersections(set<Group>& groups, 
                               set<Group>& newGroups)
    {
        for (set<Group>::iterator it2(groups.begin()); it2 != groups.end(); ++it2)
            for (set<Group>::iterator it1(groups.begin()); it1 != it2; ++it1)
            {
                DoHandleIntersections(it1, it2, newGroups);
            }
    }

    void DoHandleIntersections(set<Group>& groups1, set<Group>& groups2, 
                               set<Group>& newGroups)
    {
        for (set<Group>::iterator it2(groups2.begin()); it2 != groups2.end(); ++it2)
            for (set<Group>::iterator it1(groups1.begin()); it1 != groups1.end(); ++it1)
            {
                DoHandleIntersections(it1, it2, newGroups);
            }
    }


    bool Handle(set<Group>& oldGroups, set<Group>& newGroups)
    {
        if (oldGroups.empty() || newGroups.empty())
            return false;

        for (set<Group>::iterator itNew = newGroups.begin(); itNew != newGroups.end(); )
        {
            set<Group>::iterator itOld = oldGroups.find(*itNew);
            if (itOld != oldGroups.end())
            {
                if (itOld->countMin >= itNew->countMin && itOld->countMax <= itNew->countMax)
                {
                    newGroups.erase(itNew++);
                    continue;
                }
                else
                {
                    itNew->countMin = max(itNew->countMin, itOld->countMin);
                    itNew->countMax = min(itNew->countMax, itOld->countMax);
                    oldGroups.erase(itOld);
                }
            }
            ++itNew;
        }
        return !newGroups.empty();
    }

    void HandleIntersections(vector<Group>& groups)
    {

        set<Group> groupsSet;
        for (vector<Group>::iterator it(groups.begin()); it != groups.end(); ++it)
            Update(groupsSet, *it);
        set<Group> newGroups;

        DoHandleIntersections(groupsSet, newGroups);

        while (Handle(groupsSet, newGroups))
        {
            set<Group> newestGroups;
            DoHandleIntersections(groupsSet, newGroups, newestGroups);
            DoHandleIntersections(newGroups, newestGroups);
            groupsSet.insert(newGroups.begin(), newGroups.end());
            newGroups.swap(newestGroups);
        }

        groups.clear();
        groups.insert(groups.begin(), groupsSet.begin(), groupsSet.end());
    }

    void ProduceResult(vector<Group>& groups, set<Cell>& result)
    {
        if (!HandleUnambiguous(groups, result))
            return;

        HandleIntersections(groups);

        if (!HandleUnambiguous(groups, result))
            return;

        if (groups.empty())
        {
            results.push_back(result);
            /*
            for (int i = 0; i < 10; ++i)
            {
                for (int j = 0; j < 10; ++j)
                {
                    cout << ((result.find(Cell(i, j)) != result.end()) ? 'M' : squares[i][j]);
                }
                cout << '\n';
            }
            cout << '\n';
            */
        }
        else
        {
            vector<Group>::reverse_iterator it(groups.rbegin());
            for (int mineCount = it->countMin; mineCount <= it->countMax; ++mineCount)
            {
                if (0 == mineCount)
                {
                    vector<Group> tempGroups(groups);
                    if (Reduce(tempGroups, it->neigbors, true))
                    {
                        set<Cell> tempResult(result);
                        ProduceResult(tempGroups, tempResult);
                    }
                }
                else if (it->neigbors.size() == mineCount)
                {
                    vector<Group> tempGroups(groups);
                    if (Reduce(tempGroups, it->neigbors, false))
                    {
                        set<Cell> tempResult(result);
                        tempResult.insert(it->neigbors.begin(), it->neigbors.end());
                        ProduceResult(tempGroups, tempResult);
                    }
                }
                else
                {
                    const int nCount = it->neigbors.size();
                    int positions[8] = {0};
            	    for (int i = 0; i < mineCount; ++i)
            		    positions[i] = i;

                    for (;;)
                    {
                        vector<Cell> buffer(it->neigbors.begin(), it->neigbors.end());
                        for (int i = 0; i < mineCount; ++i)
			                if (i != positions[i])
				                swap(buffer[i], buffer[positions[i]]);

                        // stuffing
                        vector<Group> tempGroups(groups);
                        if (Reduce(tempGroups, buffer.begin(), buffer.begin() + mineCount, false)
                            && Reduce(tempGroups, buffer.begin() + mineCount, buffer.begin() + nCount, true))
                        {
                            set<Cell> tempResult(result);
                            tempResult.insert(buffer.begin(), buffer.begin() + mineCount);
                            ProduceResult(tempGroups, tempResult);
                        }

                        // next combination
		                if (positions[mineCount - 1] < nCount - 1)
			                ++positions[mineCount - 1];
		                else
		                {
                            int i;
			                for (i = mineCount - 1; --i >= 0; )
				                if (positions[i] < positions[i + 1] - 1)
					                break;
			                if (i < 0)
				                break;
			                int newValue = positions[i];
			                do
				                positions[i] = ++newValue;
			                while (++i < mineCount);
		                }
                    }
                }
            }
        }
    }
};

int _tmain(int argc, _TCHAR* argv[])
{
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER start, stop;

	QueryPerformanceCounter(&start);
    
	int runsCount = 0;

	for (;;) 
	{
        vector<Group> groups;

        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 10; ++j)
            {
                char c = squares[i][j];
                if (c >= '0' && c <= '8')
                {
                    Group group;
                    group.countMin = group.countMax = c - '0';
                    for (int i1 = max(i-1, 0); i1 <= min(i+1, 9); ++i1)
                        for (int j1 = max(j-1, 0); j1 <= min(j+1, 9); ++j1)
                            if ((i != i1 || j != j1) && (squares[i1][j1] < '0' || squares[i1][j1] > '8'))
                                group.neigbors.insert(Cell(i1, j1));

                    groups.push_back(group);
                }
            }
	//for (;;) 
	//{

        Prestidigitator prestidigitator;
        prestidigitator.ProduceResults(groups);

		runsCount++;

		QueryPerformanceCounter(&stop);
		if ((stop.QuadPart - start.QuadPart) / frequency.QuadPart < 5)
			continue;

        for (vector<set<Cell> >::iterator it(prestidigitator.results.begin())
            ; it != prestidigitator.results.end()
            ; ++it)
        {
            for (int i = 0; i < 10; ++i)
            {
                for (int j = 0; j < 10; ++j)
                {
                    cout << ((it->find(Cell(i, j)) != it->end()) ? 'M' : squares[i][j]);
                }
                cout << '\n';
            }
            cout << '\n';
        }

        //Time to solve (once): 0.0476 seconds.
		cout << "\nTime to solve (once): " <<
			double(stop.QuadPart - start.QuadPart) / frequency.QuadPart / runsCount <<
			" seconds" << endl;

		break;
	} 
	return 0;
}

