// MineSweeperRedux.cpp : Defines the entry point for the console application.
//
//	Name: Aliaksei Sanko
//	Blog: aliakseis.livejournal.com
//	Country: Belarus


#include <iostream>
#include <fstream>

#include <set>
#include <vector>

#include <algorithm>
#include <numeric>

#include <string>

#include <assert.h>
#include <time.h>

using std::ifstream;
using std::cout;
using std::endl;

using std::pair;
using std::set;
using std::vector;
using std::min;
using std::max;

using std::sort;

using std::swap;

//////////////////////////////////////////////////////////////////////////

struct Plex     // Similar to MFC CPlex
// warning variable length structure
{
	Plex* pNext;
	int dwReserved[1];    // align on 8 byte boundary

	void* data() { return this+1; }

	// like 'calloc' but no zero fill
	// may throw memory exceptions
	static Plex* Create(Plex*& pHead, size_t nMax, size_t cbElement)
	{
		assert(nMax > 0 && cbElement > 0);
		Plex* p = (Plex*) new char[sizeof(Plex) + nMax * cbElement];
				// may throw exception
		p->pNext = pHead;
		pHead = p;  // change head (adds in reverse order for simplicity)
		return p;
	}

	void FreeDataChain()       // free this one and links
	{
		Plex* p = this;
		while (p != 0)
		{
			char* bytes = (char*) p;
			p = p->pNext;
			delete[] bytes;
		}
	}
};

class AllocStack
{
	enum { BLOCK_SIZE = 10 };

	struct Node
	{
		Node* pNext;
	};

public:
	AllocStack()  
	{ 
		m_pFreeList = 0;
		m_pBlocks = 0;
	}
	~AllocStack() { RemoveAll(); }

	char* get(size_t nSize) const;
	void put(void* p) const;
	void RemoveAll()
	{
		m_pFreeList = 0;
		m_pBlocks->FreeDataChain();
		m_pBlocks = 0;
	}

private:
	mutable Node* m_pFreeList;
	mutable Plex* m_pBlocks;
};

char* AllocStack::get(size_t nSize) const
{
	if (m_pFreeList == 0)
	{
		// add another block
		Plex* newBlock = Plex::Create(m_pBlocks, BLOCK_SIZE, nSize);
		// chain them into free list
		char* pBuf = (char*) newBlock->data();
		// free in reverse order to make it easier to debug
		pBuf += nSize * (BLOCK_SIZE - 1);
		for (int i = BLOCK_SIZE - 1; i >= 0; i--, pBuf -= nSize)
		{
			Node* pNode = reinterpret_cast<Node*>(pBuf);
			pNode->pNext = m_pFreeList;
			m_pFreeList = pNode;
		}
	}
	assert(m_pFreeList != 0);  // we must have something

	char* result = reinterpret_cast<char*>(m_pFreeList);
	m_pFreeList = m_pFreeList->pNext;
	return result;
}

void AllocStack::put(void* p) const
{
	Node* pBuf = reinterpret_cast<Node*>(p);
	pBuf->pNext = m_pFreeList;
	m_pFreeList = pBuf;
}

//////////////////////////////////////////////////////////////////////////

template< typename T > struct cached_allocator : public std::allocator<T>
{
    typedef std::allocator<T> base;

    typedef T          value_type;
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;

    typedef T*         pointer;
    typedef const T*   const_pointer;

    typedef T&         reference;
    typedef const T&   const_reference;

    cached_allocator() throw() {}
    cached_allocator(const cached_allocator& a) throw() : base(a) {}
    template <typename X> cached_allocator(const cached_allocator<X>&) throw() {}
    ~cached_allocator() throw() {}

    template <typename X> struct rebind { typedef cached_allocator<X> other; };
    void deallocate(pointer p, size_type n)
    {
        if (n != 1)
            return base::deallocate(p, n);
        m_cache.put(p);
    }

    pointer allocate(size_type sz)
    {
        if (sz == 1)
        {
            return (T*)m_cache.get(sizeof T);
        }
        return base::allocate(sz);
    }

private:
    static AllocStack m_cache;
};

template <typename T>
AllocStack cached_allocator<T>::m_cache;

//////////////////////////////////////////////////////////////////////////


unsigned int FastBSF(unsigned int v)
{
    unsigned int c;     // c will be the number of zero bits on the right,
                        // so if v is 1101000 (base 2), then c will be 3
                        // NOTE: if 0 == v, then c = 31.
    if (v & 0x1)
    {
        // special case for odd v (assumed to happen half of the time)
        c = 0;
    }
    else
    {
        c = 1;
        if ((v & 0xffff) == 0)
        {
            v >>= 16;
            c += 16;
        }
        if ((v & 0xff) == 0)
        {
            v >>= 8;
            c += 8;
        }
        if ((v & 0xf) == 0)
        {
            v >>= 4;
            c += 4;
        }
        if ((v & 0x3) == 0)
        {
            v >>= 2;
            c += 2;
        }
        c -= v & 0x1;
    }
    return c;
}

typedef int Cell;

inline int ToCell(int i, int j)
{
    return i * 10 + j;
}

class CellSet
{
public:

	class iterator
	{
    public:
        iterator() : m_offset(-1), m_data(0) {}
        iterator(int offset, int data) : m_offset(offset), m_data(data) 
        {
			operator ++();
        }

        typedef std::forward_iterator_tag iterator_category;
        typedef Cell value_type;
        typedef ptrdiff_t difference_type;
        typedef Cell* pointer;
        typedef Cell& reference;

		void operator ++()
        {
            if (m_data != 0)
            {
			    m_value = m_offset + FastBSF(m_data);
			    m_data &= m_data - 1;
            }
            else
                m_offset = -1;
        }
        Cell operator *()
		{
			return m_value;
		}

		bool operator != (const iterator& other) const
        {
            return m_offset != other.m_offset || m_data != other.m_data;
        }
        bool operator == (const iterator& other) const
        {
            return m_offset == other.m_offset && m_data == other.m_data;
        }

    private:
        int m_offset; 
        unsigned int m_data;
        int m_value;
    };

    CellSet() : m_offset(0), m_data(0), m_size(0) {}

	int size() const
    {
		if (-1 == m_size)
		{
			m_size = 0;
			unsigned int v = m_data;
			if (v != 0)
			{
				v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
				v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
				m_size = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
			}
		}
        assert (m_size >= 0);
		return m_size;
    }
    bool operator < (const CellSet& other) const
    {
        return m_offset < other.m_offset 
            || m_offset == other.m_offset && m_data < other.m_data;
        //return m_data < other.m_data 
        //    || m_data == other.m_data && m_offset < other.m_offset;
    }

	iterator begin()    { return iterator(m_offset, m_data); }
    iterator end()      { return iterator(); }

    void insert(Cell value)
    {
        if (0 == m_data)
        {
            m_offset = value; 
            m_data = 1; 
            m_size = 1;
        }
        else
        {
            if (value < m_offset)
            {
                m_data <<= (m_offset - value);
                m_offset = value;
            }

            m_data |= 1 << (value - m_offset);
            m_size = -1;
        }
    }
	template<class _Iter>
	void insert(_Iter _First, _Iter _Last)
	{	// insert [_First, _Last) one at a time
	    for (; _First != _Last; ++_First)
		    insert(*_First);
	}

	size_t erase(Cell value)
    {
        if (value < m_offset || value >= m_offset + 32)
            return 0;

        int flag = 1 << (value - m_offset);
        bool found = !!(m_data & flag);
        if (found)
        {
            m_data &= ~flag;
            if (0 == m_data)
                m_size = 0;
            else
            {
                if (value == m_offset)
                {
                    int bsf = FastBSF(m_data);
                    assert(bsf != 0);
                    m_data >>= bsf;
                    m_offset += bsf;
                }
                m_size = -1;
            }
        }
        return found;
    }

    template <bool ordered>
    static bool Split(const CellSet& set1, const CellSet& set2,
        CellSet& intersect, CellSet& setI, CellSet& setJ)
    {
        if (!ordered && set2.m_offset < set1.m_offset)
        {
            return Split<true>(set2, set1, intersect, setJ, setI);
        }
        if (set2.m_offset >= set1.m_offset + 32)
            return false;

        unsigned int data2 = set1.m_data & (set2.m_data << (set2.m_offset - set1.m_offset));

        intersect = CellSet(set1.m_offset, set1.m_data & data2);
        if (intersect.size() < 2)
            return false;

        setI = CellSet(set1.m_offset, set1.m_data & ~data2);
        setJ = CellSet(set2.m_offset, set2.m_data & ~(set1.m_data >> (set2.m_offset - set1.m_offset)));
        return true;
    }

    bool operator == (const CellSet& other)
    {
        return m_data == other.m_data && (0 == m_data || m_offset == other.m_offset);
    }

    int GetOffset() const { return m_offset; }

private:
    CellSet(int offset, unsigned int data) : m_offset(offset), m_data(data), m_size(-1) 
    {
        if (m_data != 0)
        {
            int shift = FastBSF(m_data);
            m_data >>= shift;
            m_offset += shift;
        }
        else
            m_size = 0;
    }

    int m_offset;
    unsigned int m_data;
	mutable int m_size;
};

struct Group
{
    CellSet neigbors;
    int countMin, countMax;

    bool operator < (const Group& other) const
    {
        return neigbors < other.neigbors;
    }
};


typedef std::set<Group, std::less<Group>, cached_allocator<Group> > GroupSet;


class Prestidigitator
{
public:
    vector<set<Cell> > results;

    void ProduceResults(const vector<Group>& groups)
    {
        vector<Group> temp(groups);
        set<Cell> result;
        ProduceResult(temp, result);
    }

private:
    template<typename T>
    bool Reduce(vector<Group>& groups, T begin, T end, const bool areZeros)
    {
        for (int i = (int)groups.size(); --i >= 0;)
        {
            Group& group = groups[i];
            for (T it = begin; it != end; ++it)
                if (group.neigbors.erase(*it))
                    if (!areZeros)
                    {
                        if (group.countMin > 0)
                            group.countMin--;
                        group.countMax--;
                    }
                    else
                    {
                        if (group.countMax > (int)group.neigbors.size())
                            group.countMax = (int)group.neigbors.size();
                    }

            if (group.countMin > group.countMax)
                return false;

            if (0 == group.neigbors.size())
                groups.erase(groups.begin() + i);
        }

        return true;
    }

    template<typename T>
    bool Reduce(vector<Group>& groups, T& cells, const bool areZeros)
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
                    CellSet cells(it->neigbors);
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

    bool Update(GroupSet& groups, const Group& newGroup)
    {
        pair<GroupSet::iterator, bool> inserted = groups.insert(newGroup);
        if (inserted.second)
            return true;

        bool result = false;
        if (inserted.first->countMin < newGroup.countMin)
        {
            result = true;
            const_cast<Group&>(*inserted.first).countMin = newGroup.countMin;
        }
        if (inserted.first->countMax > newGroup.countMax)
        {
            result = true;
            const_cast<Group&>(*inserted.first).countMax = newGroup.countMax;
        }

        return result;
    }


    void DoHandleIntersections(GroupSet::iterator it1, GroupSet::iterator it2, GroupSet& newGroups)
    {
        CellSet intersect, setI, setJ;
        if (!CellSet::Split<false>(it1->neigbors, it2->neigbors, intersect, setI, setJ))
            return;

        int intersectMax = min((int)intersect.size(), min(it1->countMax, it2->countMax));

        int intersectMin = max(0, 
            (int)(intersect.size() - min(it1->neigbors.size() - it1->countMin, it2->neigbors.size() - it2->countMin)));

        if (setI.size() > 1)
        {
            Group groupI;
            groupI.neigbors = setI;
            groupI.countMin = max(it1->countMin - intersectMax, 0);
            groupI.countMax = min(it1->countMax - intersectMin, (int)groupI.neigbors.size());
            Update(newGroups, groupI);
        }
        if (setJ.size() > 1)
        {
            Group groupJ;
            groupJ.neigbors = setJ;
            groupJ.countMin = max(it2->countMin - intersectMax, 0);
            groupJ.countMax = min(it2->countMax - intersectMin, (int)groupJ.neigbors.size());
            Update(newGroups, groupJ);
        }

        Group group;
        group.neigbors = intersect;
        group.countMin = intersectMin;
        group.countMax = intersectMax;
        Update(newGroups, group);
    }

    void DoHandleIntersections(GroupSet& groups, GroupSet& newGroups)
    {
        GroupSet::iterator itEnd(groups.end());
        for (GroupSet::iterator it1(groups.begin()); it1 != itEnd; ++it1)
        {
            const int border = it1->neigbors.GetOffset() + 23;
            for (GroupSet::iterator it2(it1); ++it2 != itEnd; )
            {
                if (it2->neigbors.GetOffset() > border)
                    break;
                DoHandleIntersections(it1, it2, newGroups);
            }
        }
    }

    void DoHandleIntersections(GroupSet& groups1, GroupSet& groups2, GroupSet& newGroups)
    {
        for (GroupSet::iterator it1end(groups1.end()), it1(groups1.begin()); it1 != it1end; ++it1)
        {
            const int border = it1->neigbors.GetOffset() + 23;
            for (GroupSet::iterator it2end(groups2.end()), it2(groups2.begin()); it2 != it2end; ++it2)
            {
                if (it2->neigbors.GetOffset() > border)
                    break;
                DoHandleIntersections(it1, it2, newGroups);
            }
        }
    }


    bool Handle(GroupSet& oldGroups, GroupSet& newGroups)
    {
        if (oldGroups.empty() || newGroups.empty())
            return false;

        for (GroupSet::iterator itNew = newGroups.begin(); itNew != newGroups.end(); )
        {
            GroupSet::iterator itOld = oldGroups.find(*itNew);
            if (itOld != oldGroups.end())
            {
                if (itOld->countMin >= itNew->countMin && itOld->countMax <= itNew->countMax)
                {
                    newGroups.erase(itNew++);
                    continue;
                }
                else
                {
                    const_cast<Group&>(*itNew).countMin = max(itNew->countMin, itOld->countMin);
                    const_cast<Group&>(*itNew).countMax = min(itNew->countMax, itOld->countMax);
                    oldGroups.erase(itOld);
                }
            }
            ++itNew;
        }
        return !newGroups.empty();
    }

    void HandleIntersections(vector<Group>& groups)
    {
        GroupSet groupsSet;
        for (vector<Group>::iterator it(groups.begin()); it != groups.end(); ++it)
            Update(groupsSet, *it);
        GroupSet newGroups;

        DoHandleIntersections(groupsSet, newGroups);

        while (Handle(groupsSet, newGroups))
        {
            GroupSet newestGroups;
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

        if (!groups.empty())
        {
            HandleIntersections(groups);

            if (!HandleUnambiguous(groups, result))
                return;
        }

        if (groups.empty())
        {
            results.push_back(result);
        }
        else
        {
            vector<Group>::iterator it(groups.begin());
            //vector<Group>::reverse_iterator it(groups.rbegin());
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

int main(int /*argc*/, char* argv[])
{
	// The input file should be in the same location as our executable. 
    std::string path("input.txt");
    {
        std::string dir(argv[0]);
        std::string::size_type pos = dir.find_last_of("\\/");
        if (std::string::npos != pos)
        {
            path = dir.substr(0, pos + 1) + path;
        }
    }

	ifstream inFile(path);

    char squares[10][11] = {0};
    for (int i = 0; i < 10; ++i)
        inFile.getline(squares[i], 11);

    clock_t start = clock();
    
	int runsCount = 0;
    
    //*
	for (;;) 
	{
    //*/
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
                                group.neigbors.insert(ToCell(i1, j1));

                    groups.push_back(group);
                }
            }
    /*
	for (;;) 
	{
    //*/
        Prestidigitator prestidigitator;
        prestidigitator.ProduceResults(groups);

		runsCount++;

        clock_t stop = clock();

		if ((stop - start) / CLOCKS_PER_SEC < 5)
			continue;

        for (vector<set<Cell> >::iterator it(prestidigitator.results.begin())
            ; it != prestidigitator.results.end()
            ; ++it)
        {
            for (int i = 0; i < 10; ++i)
            {
                for (int j = 0; j < 10; ++j)
                {
                    cout << ((it->find(ToCell(i, j)) != it->end()) ? 'M' : squares[i][j]);
                }
                cout << '\n';
            }
            cout << '\n';
        }

        //Time to solve (once): 0.0476 seconds.
		cout << "Time to solve (once): " <<
			double(stop - start) / CLOCKS_PER_SEC / runsCount <<
			" seconds" << endl;

		break;
	} 
	return 0;
}
