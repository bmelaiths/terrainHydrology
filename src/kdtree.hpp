#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <queue>

#include "point.hpp"

template <typename T>
struct Node
{
    Node *left = NULL, *right = NULL;
    Point loc;
    T idx;
    //these two numbers are used for the lock hierarchy
    //size_t idOfCreatorThread, id;
    //there should also be some kind of lock

    Node(Point loc, T idx)
    : loc(loc), idx(idx)
    {};
};

template <typename T>
class KDTree
{
    private:
    Node<T>* root = NULL;
    std::vector<Node<T>*> allNodes; // this will be good for the destructor

    private:
    Node<T>* tryInsert(Node<T>* finalNode, bool isXlevel, Node<T>* toIns);

    public:
    void insert(Point loc, T idx);
    std::vector<T> rangeSearch(Point loc, float radius);
    std::vector<T> breadthFirstSearch();
};

template <typename T>
Node<T>* KDTree<T>::tryInsert(Node<T>* finalNode, bool isXlevel, Node<T>* toIns)
{
    //if this is an X level
        //if loc.x is less than finalNode's x
            //if there is no left child, then this IS final level
            //else NOT
        //else
            //if there is no right child, then this IS final level
            //else NOT

    //if this is a Y level
        //if loc.y is less than finalNode's y
            //if there is no left child, then this IS final level
            //else NOT
        //else
            //if there is no right child, then this IS final level
            //else NOT
    
    //if this is an x level, then we compare x values, else y values
    float nodeCompValue = isXlevel ? finalNode->loc.x() : finalNode->loc.y();
    float insCompValue  = isXlevel ?     toIns->loc.x() : toIns->loc.y()    ;

    if (insCompValue < nodeCompValue)
    { //work on the left side of the tree
        if (finalNode->left == NULL)
        { //if the new node can be inserted, do so
            finalNode->left = toIns;
        }
        return finalNode->left;
    }
    else
    { //work on the right side of the tree
        if (finalNode->right == NULL)
        {
            finalNode->right = toIns;
        }
        return finalNode->right;
    }
}

template <typename T>
void KDTree<T>::insert(Point loc, T idx)
{
    Node<T>* newNode = new Node<T>(loc, idx);

    if (root == NULL)
    {
        root = newNode;
        return;
    }
    

    Node<T>* finalNode = root;
    Node<T>* nextNode = NULL;
    bool isXlevel = true;
    //tryInsert() will return a pointer to the node at the next level. If there
    //wasn't a node there before, then it will insert newNode and return a
    //pointer to it, since it is the node at the next level. This loop will run
    //until that happens
    while ((nextNode=tryInsert(finalNode, isXlevel, newNode)) != newNode)
    {
        finalNode = nextNode;
        isXlevel = !isXlevel;
    }
    //loop invariant: finalNode points to the node that the new node will be
    //inserted under (unless the next iteration finds another level)

    allNodes.push_back(newNode);
}

template <typename T>
struct rangeSearchStruct
{
    Node<T> *node;
    bool isXlevel;
};
template <typename T>
std::vector<T> KDTree<T>::rangeSearch(Point loc, float radius)
{
    std::vector<T> foundIdxes;

    std::queue<rangeSearchStruct<T>> toSearch;
    //insert the root node, and mark it as a level sorted by X
    toSearch.push({.node=root, .isXlevel=true});
    while (toSearch.size() > 0)
    {
        //extract a node to evaluate
        Node<T> *v    = toSearch.front().node;
        bool isXlevel = toSearch.front().isXlevel;
        toSearch.pop();

        //if the node doesn't exist, don't evaluate it
        if (v == NULL) {
            continue;
        }

        //figure out which values to compate
        float nodeCompValue, otherNodeCompValue;
        float queryCompValue, otherQueryCompValue;
        if (isXlevel)
        { //if this is an X level, evaluate the x values
            nodeCompValue = v->loc.x();
            otherNodeCompValue = v->loc.y();
            queryCompValue = loc.x();
            otherQueryCompValue = loc.y();
        }
        else
        {
            nodeCompValue = v->loc.y();
            otherNodeCompValue = v->loc.x();
            queryCompValue = loc.y();
            otherQueryCompValue = loc.x();
        }

        //based on the relevant dimension, decide what to do
        if (nodeCompValue < (queryCompValue - radius))
        { //the node cannot be in range, but some of its right children could be
            toSearch.push({.node=v->right, .isXlevel=!isXlevel});
        }
        else if (
            (queryCompValue - radius) <= nodeCompValue &&
            nodeCompValue < (queryCompValue + radius))
        {
            //if the relevant dimension is in range, its irrelevant dimension could be, too
            if
            (
                (otherQueryCompValue - radius) <= otherNodeCompValue &&
                 otherNodeCompValue < (otherQueryCompValue + radius)
            )
            {
                foundIdxes.push_back(v->idx);
            }
            //descendants in both trees should be searched
            toSearch.push({.node=v->left , .isXlevel=!isXlevel});
            toSearch.push({.node=v->right, .isXlevel=!isXlevel});
        }
        else
        {
            toSearch.push({.node=v->left, .isXlevel=!isXlevel});
        }
    }
    
    return foundIdxes;
}

template <typename T>
std::vector<T> KDTree<T>::breadthFirstSearch() {
    std::vector<T> discovered;
    
    std::queue<Node<T>*> q;

    discovered.push_back(root->idx);
    q.push(root);
    while (q.size() > 0)
    {
        Node<T>* v = q.front();
        q.pop();
        if (v->left != NULL) {
            discovered.push_back(v->left->idx);
            q.push(v->left);
        }
        if (v->right != NULL) {
            discovered.push_back(v->right->idx);
            q.push(v->right);
        }
    }

    return discovered;
}

#endif