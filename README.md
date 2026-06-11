![Image](https://github.com/user-attachments/assets/7f96012f-c05b-4b85-8d43-531dff6bed3a)

# MineSweeperRedux

A high-performance C++ Minesweeper constraint solver that enumerates all valid mine configurations for a partially revealed board.

This project demonstrates advanced constraint propagation, set intersection analysis, combinatorial search, and custom memory allocation techniques to efficiently solve Minesweeper puzzles.

## Overview

MineSweeperRedux accepts a partially solved 10×10 Minesweeper board and computes every possible mine placement consistent with the revealed numeric clues.

Unlike simple Minesweeper assistants that identify only immediately safe moves, this solver performs a complete logical and combinatorial analysis of the puzzle state, generating all valid solutions.

The project was originally developed as an experiment in:

* Constraint satisfaction problems (CSP)
* Efficient set representation
* Search-space reduction
* Custom STL allocators
* Memory pool optimization
* Bitwise algorithms

---

## Features

### Constraint-Based Solving

Each revealed Minesweeper number is converted into a constraint:

> Among these neighboring unrevealed cells, exactly N mines exist.

The solver builds a system of overlapping constraints and repeatedly refines them until all possibilities are resolved.

### Automatic Deduction

The engine immediately detects:

* Guaranteed mines
* Guaranteed safe cells
* Contradictory states
* Fully resolved groups

Examples:

```
{A,B,C} = 0
```

means all cells are safe.

```
{A,B,C} = 3
```

means all cells are mines.

These deductions are propagated through all remaining constraints.

---

### Intersection Analysis

A major optimization is the analysis of overlapping constraint groups.

Example:

```
{A,B,C,D} = 2
{C,D,E,F} = 2
```

The solver computes:

```
Intersection: {C,D}
Left:         {A,B}
Right:        {E,F}
```

and derives tighter bounds for each subset.

This frequently eliminates large portions of the search tree before any brute-force enumeration begins.

---

### Complete Solution Enumeration

When pure deduction is no longer sufficient, the solver recursively explores all remaining valid combinations.

Every generated solution satisfies all board constraints.

Output consists of all valid mine layouts.

---

## Internal Architecture

### CellSet

`CellSet` is a compact bitset-like container optimized for Minesweeper neighborhoods.

Characteristics:

* Stores up to 32 cells
* Uses a single 32-bit integer
* Constant-time insertion/removal
* Fast iteration using bit scanning
* Extremely low memory overhead

Example:

```cpp
CellSet neighbors;
neighbors.insert(cell);
```

---

### FastBSF

The project includes a custom implementation of:

```cpp
unsigned int FastBSF(unsigned int v);
```

which performs a fast "Bit Scan Forward" operation.

Used for:

* Iteration over active bits
* Compact set traversal
* Efficient cell extraction

This avoids heavier generic containers during intensive search operations.

---

### Group

A `Group` represents a Minesweeper constraint.

```cpp
struct Group
{
    CellSet neighbors;
    int countMin;
    int countMax;
};
```

Meaning:

```
countMin ≤ mines(neighbors) ≤ countMax
```

Initially:

```
countMin == countMax
```

but intersection processing may relax or tighten these bounds.

---

### Prestidigitator

The core solving engine.

Responsibilities:

1. Constraint reduction
2. Propagation of known mines/safe cells
3. Constraint intersection analysis
4. Recursive search
5. Solution generation

Main entry point:

```cpp
Prestidigitator::ProduceResults(...)
```

which produces all valid mine configurations.

---

## Memory Optimization

### Cached Allocator

The solver uses a custom STL allocator:

```cpp
cached_allocator<T>
```

designed to accelerate frequent allocation of small objects stored inside STL containers.

Benefits:

* Reduced heap fragmentation
* Lower allocation overhead
* Improved search performance

---

### AllocStack

A lightweight object pool used by the allocator.

Features:

* Block allocation
* Free-list recycling
* Constant-time allocation/deallocation

This is especially useful because constraint generation creates many temporary objects.

---

### Plex

Inspired by MFC's `CPlex`.

Provides:

* Chunked memory allocation
* Linked allocation blocks
* Fast bulk cleanup

Instead of individually freeing thousands of small objects, entire memory chains can be released at once.

---

## Input Format

The solver expects a file named:

```text
input.txt
```

located next to the executable.

Example:

```text
1.........
2.........
..........
..........
..........
..........
..........
..........
..........
..........
```

### Cell Meaning

| Character           | Meaning       |
| ------------------- | ------------- |
| 0-8                 | Revealed clue |
| Any other character | Unknown cell  |

Unknown cells are treated as potential mine locations.

---

## Output

For every valid solution:

```text
1MM......
2M.......
.........
...
```

where:

```text
M
```

marks a mine.

After enumeration, timing statistics are displayed:

```text
Time to solve (once): 0.0476 seconds
```

---

## Algorithm Summary

1. Parse board.
2. Generate constraints from revealed numbers.
3. Apply deterministic reductions.
4. Analyze overlapping constraints.
5. Repeat until no new information appears.
6. Recursively enumerate remaining possibilities.
7. Output every valid solution.

---

## Complexity

Worst-case complexity remains exponential, as Minesweeper solving is closely related to NP-complete constraint satisfaction problems.

However, several optimizations dramatically reduce practical runtime:

* Compact bitset representation
* Constraint propagation
* Group intersection deduction
* Search-space pruning
* Custom allocator
* Memory pooling

For realistic boards, the solver is significantly faster than naïve brute-force enumeration.

---

## Building

Any modern C++ compiler supporting STL containers should build the project.

Example:

```bash
g++ MineSweeperRedux.cpp -O3 -std=c++11 -o MineSweeperRedux
```

or using Visual Studio:

```bash
cl /O2 MineSweeperRedux.cpp
```

---

## Educational Value

This project is an interesting example of:

* Constraint satisfaction systems
* Logical inference engines
* Bit manipulation techniques
* Search optimization
* Custom allocators
* Memory pooling
* Efficient STL integration

It demonstrates how a seemingly simple puzzle can be transformed into a sophisticated constraint-solving problem.

