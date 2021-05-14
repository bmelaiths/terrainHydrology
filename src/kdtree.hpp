#ifndef KDTREE_H
#define KDTREE_H

#include <queue>
#include <vector>
#include <algorithm>
#include <math.h>
#include <omp.h>
#include <cfloat>
#include <stdio.h>

#include "point.hpp"

/**
 * @brief A node for use in a 2DTree
 * 
 * @tparam T The type of the data to store
 */
template <typename T>
struct Node
{
    Node *left = NULL, *right = NULL;
    /**
     * @brief The location of the data
     */
    Point loc;
    /**
     * @brief The data
     */
    T idx;
    bool isXlevel;

    /**
     * @brief Construct a new Node object
     * 
     * @param loc The location of the data
     * @param idx The data stored in the point
     */
    Node(Point loc, T idx)
    : loc(loc), idx(idx)
    {}
    ~Node()
    {
        if (left != NULL)
        {
            delete left;
        }
        if (right != NULL)
        {
            delete right;
        }
    }
    Node(const Node& other)
    : left(other.left), right(other.right), loc(other.loc),
      idx(other.idx), isXlevel(other.isXlevel)
    {}
    Node(Node&& other) noexcept
    : left(std::move(other.left)), right(std::move(other.right)),
      loc(std::move(other.loc)), idx(std::move(other.idx)),
      isXlevel(std::move(other.isXlevel))
    {
        other.left = NULL;
        other.right = NULL;
    }

    Node& operator=(const Node& other)
    {
        if (this == &other)
        {
            return *this;
        }

        left = other.left;
        right = other.right;
        loc = other.loc;
        idx = other.idx;
        isXlevel = other.isXlevel;

        return *this;
    }
    Node& operator=(Node&& other) noexcept
    {
        if (this == &other)
        {
            return *this;
        }

        left = std::move(other.left);
        right = std::move(other.right);
        loc = std::move(other.loc);
        idx = std::move(other.idx);
        isXlevel = std::move(other.isXlevel);

        other.left = NULL;
        other.right = NULL;

        return *this;
    }
};

/**
 * @brief A binary tree for storing 2D points and associated data
 * 
 * @tparam T The data to store
 */
template <typename T>
class KDTree
{
    private:
    Node<T>* root = NULL;
    std::vector<Node<T>*> allNodes;

    private:
    /**
     * @brief Tries to insert the point as a child of `finalNode`
     * 
     * This method will return a pointer to the node that is at the next
     * level. If `toIns` belongs as a child of `finalNode`, then a pointer
     * to `toIns` will be returned, otherwise, a pointer to the relevant
     * child of `finalNode` will be returned, and this method can be called
     * with that node.
     * 
     * @param finalNode The node that _may_ be the rightful parent of `toIns`
     * @param isXlevel Specifies whether this level of the tree is sorted by X (if true), or Y (if false)
     * @param toIns The node that you want to insert into the tree
     * @return Node<T>* Either `toIns`, or the child of finalNode to follow
     */
    Node<T>* tryInsert(Node<T>* finalNode, bool isXlevel, Node<T>* toIns);
    /**
     * @brief Returns the node that should be the root node of a range
     * 
     * Used to reconstruct the tree. The user specifies a range within the
     * vector of all nodes, from [begin,end), and constructs a balanced
     * 2DTree that includes all of those nodes.
     * 
     * @param begin The start index of the allNodes vector
     * @param end The end index (not inclusive) of the allNodes vector
     * @param isXlayer Should the first layer of this tree sort by X value
     * @return Node<T>* Returns the root node of the new tree
     */
    Node<T>* makeRoot(size_t begin, size_t end, bool isXlayer);
    /**
     * @brief Gets the maximum depth of the tree
     * 
     * Used to determine whether or not the tree is balanced
     * 
     * @param node The node to start from
     * @return size_t The maximum depth of the tree
     */
    size_t maxDepth(Node<T>* node);
    /**
     * @brief Makes a copy of the subtree rooted at `node`
     * 
     * @param node The root of the subtree to copy
     * @return Node<T>* The root node of a copy of the subtree
     */
    Node<T>* copySubTree(Node<T> *node);

    public:
    KDTree();
    ~KDTree();
    KDTree(const KDTree<T>& other);
    KDTree(KDTree<T>&& other) noexcept;

    KDTree& operator=(const KDTree& other);
    KDTree& operator=(KDTree&& other) noexcept;

    /**
     * @brief Insert a point and its associated data into the tree
     * 
     * @param loc The location of the data
     * @param idx The data of the location
     */
    void insert(Point loc, T idx);
    /**
     * @brief Get all data within an area
     * 
     * The area searched is square, with its center `loc`, and
     * 2 * radius wide and high.
     * 
     * @param loc The center of the area to search
     * @param radius The width and height of the area to search
     * @return std::vector<T> All the data situated in that area
     */
    std::vector<T> searchRange(Point loc, float radius);
    /**
     * @brief Gets all data in the tree in BFS order
     * 
     * This is mostly used for testing and debugging
     * 
     * @return std::vector<T> All data in BFS order
     */
    std::vector<T> breadthFirstSearch();
    /**
     * @brief The maximum depth of the tree.
     * 
     * Used to determine how balanced the tree is
     * 
     * @return size_t 
     */
    size_t getDepth();
    /**
     * @brief Reconstruct the tree to be as balanced as possible
     * 
     */
    void reconstruct();
};

template <typename T>
KDTree<T>::KDTree()
: root(NULL)
{}

template <typename T>
KDTree<T>::~KDTree()
{
    if (root != NULL)
    {
        delete root;
    }
}

template <typename T>
KDTree<T>::KDTree(const KDTree<T>& other)
{
    root = copySubTree(other.root);
    allNodes = other.allNodes;
}

template <typename T>
KDTree<T>::KDTree(KDTree<T>&& other) noexcept
: root(std::move(other.root)), allNodes(std::move(allNodes))
{
    other.root = NULL;
}

template <typename T>
KDTree<T>& KDTree<T>::operator=(const KDTree<T>& other)
{
    if (this == &other)
    {
        return *this;
    }
    root = copySubTree(other.root);
    allNodes = other.allNodes;
}

template <typename T>
KDTree<T>& KDTree<T>::operator=(KDTree<T>&& other) noexcept
{
    root = std::move(other.root);
    allNodes = std::move(allNodes);

    other.root = NULL;
}

template <typename T>
Node<T>* KDTree<T>::copySubTree(Node<T> *node)
{
    // this method operates in a recursive manner
    // it makes a copy of this node (using the Node() copy constructor)
    // and calls this method for its children
    // if this node is NULL, then the method will return NULL
    if (node == NULL)
    {
        return NULL;
    }
    Node<T> *nodeCopy = new Node<T>(*node);
    nodeCopy->left = copySubTree(node->left);
    nodeCopy->right = copySubTree(node->right);
    return nodeCopy;
}

template <typename T>
size_t KDTree<T>::maxDepth(Node<T>* node)
{
    // this method operates in a recursive manner
    // it calculates the maximum depth of the tree backwards
    // the 'depth' of a leaf node is 1
    // the depth of a node with one or two leaf children is 2
    // but the depth of a node is the MAXIMUM depth of either of its children
    int max = 1;
    if (node->left != NULL)
    {
        max = 1 + maxDepth(node->left);
    }
    if (node->right != NULL)
    {
        int rightDepth = 1 + maxDepth(node->right);
        max = (rightDepth > max) ? rightDepth : max;
    }
    return max;
}

template <typename T>
size_t KDTree<T>::getDepth()
{
    return maxDepth(root);
}

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
            toIns->isXlevel = !isXlevel;
        }
        return finalNode->left;
    }
    else
    { //work on the right side of the tree
        if (finalNode->right == NULL)
        {
            finalNode->right = toIns;
            toIns->isXlevel = !isXlevel;
        }
        return finalNode->right;
    }
}

template <typename T>
void KDTree<T>::insert(Point loc, T idx)
{
    // if the tree is unbalanced, then reconstruct the tree
    if (allNodes.size() > 50)
    {
        float inefficiency =
            getDepth() / (log(allNodes.size()) / log(2))
        ;
        if (inefficiency > 2)
        {
            reconstruct();
        }
    }

    Node<T>* newNode = new Node<T>(loc, idx);
    allNodes.push_back(newNode);

    // if this tree is empty, insert the node as the root of the tree
    if (root == NULL)
    {
        newNode->isXlevel = true;
        root = newNode;
        return;
    }

    Node<T>* finalNode = root;
    Node<T>* nextNode = NULL;
    bool isXlevel = true;
    //tryInsert() will return a pointer to the node at the next level. If there
    //wasn't a node there before, then it will insert newNode and return a
    //pointer to it, since it is now the node at the next level. This loop will run
    //until that happens
    while ((nextNode=tryInsert(finalNode, isXlevel, newNode)) != newNode)
    {
        finalNode = nextNode;
        isXlevel = !isXlevel;
    }
    //loop invariant: finalNode points to the node that the new node will be
    //inserted under (unless the next iteration finds another level)
}

template <typename T>
struct rangeSearchStruct
{ // simple struct used to aid in searching
    Node<T> *node;
    bool isXlevel;
};
template <typename T>
std::vector<T> KDTree<T>::searchRange(Point loc, float radius)
{
    std::vector<T> foundIdxes;

    if (root == NULL)
    {
        return foundIdxes;
    }

    std::queue<rangeSearchStruct<T>> toSearch;
    //insert the root node, and mark it as a level sorted by X
    toSearch.push({.node=root, .isXlevel=true});
    // toSearch.push({.node=root, .isXlevel=true});
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
    
    /* Simple BFS algorithm */

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

template <typename T>
class CompareNodesX
{
    public:
    bool operator() (Node<T>* first, Node<T>* second)
    {
        return first->loc.x() < second->loc.x();
    }
};

template <typename T>
class CompareNodesY
{
    public:
    bool operator() (Node<T>* first, Node<T>* second)
    {
        return first->loc.y() < second->loc.y();
    }
};

template <typename T>
Node<T>* KDTree<T>::makeRoot(size_t begin, size_t end, bool isXlevel)
{
    /* This method will sort a subtree, as delineated by 'begin'
     * and 'end', which are indices within the 'allNodes' vector.
     */

    if (begin == end)
    { // if there are no nodes within this tree, return NULL
        return NULL;
    }
    
    if (end-begin == 1)
    { // if there is only one node within this tree, then this
      //node is the root, and it has no children
        allNodes[begin]->left = NULL;
        allNodes[begin]->right = NULL;
        allNodes[begin]->isXlevel = isXlevel; // this is SUPER DUPER POOPER SCOOPER important
        return allNodes[begin];
    }
    
    //sort this subtree by the appropriate dimension
    if (isXlevel)
    {
        std::sort(
            allNodes.begin() + begin,
            allNodes.begin() + end,
            CompareNodesX<T>()
        );
    }
    else
    {
        std::sort(
            allNodes.begin() + begin,
            allNodes.begin() + end,
            CompareNodesY<T>()
        );
    }
    
    //the root of this subtree is the median of nodes, along the appropriate dimension
    Node<T>* newRoot =  allNodes[begin + (int)(end-begin)/2];
    newRoot->isXlevel = isXlevel; // this is SUPER important

    //the child nodes should be adequately-constructed subtrees
    newRoot->left =  makeRoot(begin, begin + (end-begin)/2, !isXlevel);
    newRoot->right = makeRoot((int) begin + (end-begin)/2 + 1, end, !isXlevel);

    return newRoot;
}

template <typename T>
void KDTree<T>::reconstruct()
{
    //sort all nodes by the appropriate dimension, which, for the root
    //level, is X
    std::sort(allNodes.begin(), allNodes.end(), CompareNodesX<T>());

    //the root of this tree should be the median of nodes
    root = allNodes[(int)allNodes.size()/2];
    root->isXlevel = true; // this is also SUPER important

    //the two subtrees should be adequately constructed subtrees
    root->left = makeRoot(0, (int) allNodes.size()/2, false);
    root->right = makeRoot((int) allNodes.size()/2 + 1, allNodes.size(), false);
}

#endif