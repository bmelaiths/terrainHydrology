#ifndef KDTREE_H
#define KDTREE_H

#include <vector>

#include "point.hpp"

struct Node
{
    Node *left = NULL, *right = NULL;
    Point loc;
    size_t idx;
    //these two numbers are used for the lock hierarchy
    //size_t idOfCreatorThread, id;
    //there should also be some kind of lock

    Node(Point loc, size_t idx)
    : loc(loc), idx(idx)
    {};
};

class KDTree
{
    private:
    Node* root;
    std::vector<Node*> allNodes; // this will be good for the destructor

    private:
    Node* tryInsert(Node* finalNode, bool isXlevel, Node* toIns);

    public:
    KDTree(std::vector<Point> initPoints, std::vector<size_t> initIdxes);
    ~KDTree();

    public:
    void insert(Point loc, size_t idx);
    std::vector<size_t> rangeSearch(Point loc, float radius);
    std::vector<size_t> breadthFirstSearch();
};


#endif