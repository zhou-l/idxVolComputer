// From the software distribution accompanying the textbook
// "A Practical Introduction to Data Structures and Algorithm Analysis,
// Third Edition (C++)" by Clifford A. Shaffer.
// Source code Copyright (C) 2007-2011 by Clifford A. Shaffer.
#ifndef KDTREE_H
#define KDTREE_H
#include "global.h"
#include "BSTNode.h"
#include "book.h"
#include <fstream>
using std::ofstream;
using std::vector;
// k-d tree class
template <typename E> 
class KD {
public:
	// Constructor
	KD(int Din) { root = NULL; D = Din; nodecount = 0; }

	void insert(const E& it) {
		root = inserthelp(root, it, 0);
		nodecount++;
	}

	void clear()
	{
		clearhelp(root); root = NULL; nodecount = 0;
	}

	void remove(int* coord, E& it) {
		BinNode<E>* t = NULL;
		root = removehelp(root, coord, t, 0);
		if (t != NULL) {
			it = t->element();
			nodecount--;
			delete t;
		}
	}

	bool find(E& e, int* coord) const
	{
		return findhelp(root, coord, e, 0);
	}

	void regionsearch(int* coord, int radius) const {
		regionhelp(root, coord, radius, 0);
	}
	// Rectangular queries
	void queryRegion(int* startCoord, int* endCoord) const{
		queryRegionHelp(root, startCoord, endCoord, 0);
	}
	
	void queryRegion(float* startCoord, float* endCoord) const{
		cout << "Use float!" << endl;
		queryRegionHelp(root, startCoord, endCoord, 0);
	}
	// return query result
	void queryRegionGetResult(float* startCoord, float* endCoord, std::vector<E>& result)
	{
		cout << "Query with result as output!" << endl;
		queryRegionGetResultHelp(root, startCoord, endCoord, 0, result);
	}
	// Balance tree construction
	void buildBalanceTree(std::vector<E> list)
	{
		root = balanceBuildHelp(root, list, 0);
	}

	int size() { return nodecount; }

	void print() const {
		if (root == NULL) cout << "The kd-tree is empty.\n";
		else {
			printhelp(root, 0);
			cout << endl;
		}
	}
	
	void printToFile(ofstream& ofFile) const{
		if (root == NULL)
			ofFile << "The kd-tree is empty.\n";
		else{
			printToFileHelp(root, ofFile, 0);
			ofFile << endl;
		}

	}

private:
	BinNode<E>* root;
	int D;                 // Dimension of the key
	int nodecount;         // Number of nodes in the KD tree

	// Private "helper" functions

	// Return true iff the two points are within radius distance of each other
	bool InCircle(int* c1, int* c2, int radius) const {
		int total = 0;
		for (int i = 0; i<D; i++) // Compute sum of squares
			total += (c1[i] - c2[i]) * (c1[i] - c2[i]);
		return (total <= (radius*radius));
	}

	// return true iff the point is inside/on the rectangle defined by [startCoord, endCoord]
	bool InRect(int* c, int* startCoord, int* endCoord) const {
		bool insideRect = true;
		for (int i = 0; i < D; i++)
		{
			if (c[i] < startCoord[i] || c[i] > endCoord[i])
			{
				insideRect = false;
				break;
			}
		}
		return insideRect;
	}

	bool InRect(float* c, float* startCoord, float* endCoord) const {
		bool insideRect = true;
		for (int i = 0; i < D; i++)
		{
			if (c[i] < startCoord[i] || c[i] > endCoord[i])
			{
				insideRect = false;
				break;
			}
		}
		return insideRect;
	}

	// True iff two coords are equal
	bool EqualCoord(const int* a, const int* b) const {
		for (int i = 0; i<D; i++)
		if (a[i] != b[i]) return false;
		return true;
	}

	// True iff two coords are equal
	bool EqualCoord(const float* a, const float* b) const {
		for (int i = 0; i<D; i++)
		if (a[i] != b[i]) return false;
		return true;
	}

	// Return the node with the smaller value for this discrim
	BinNode<E>* min(BinNode<E>* node1,
		BinNode<E>* node2, int discrim) const {
		if (node1 == NULL) return node2;
		if (node2 == NULL) return node1;
		if (((node2->element())->coord())[discrim] <
			((node1->element())->coord())[discrim])
			return node2;
		else return node1;
	}

	// Clean up the KD tree by releasing space back free store
	void clearhelp(BinNode<E>* root) {
		if (root == NULL) return;
		clearhelp(root->left());
		clearhelp(root->right());
		delete root;
	}

	// Insert a node into the KD tree, returning the updated tree
	BinNode<E>* inserthelp(BinNode<E>* root,
		const E& it, int discrim) {
		if (root == NULL)  // Empty tree: create node
			return (new BSTNode<int, E>(0, it, NULL, NULL));
		if (((root->element())->coord())[discrim] > (it->coord())[discrim])
			root->setLeft(inserthelp(root->left(), it, (discrim + 1) % D));
		else root->setRight(inserthelp(root->right(), it, (discrim + 1) % D));
		return root;       // Return tree with node inserted
	}

	// Create balance KD-tree
	BinNode<E>* balanceBuildHelp(BinNode<E>* root, vector<E>& list, int discrim)
	{
		size_t id;
		if (list.empty())
			return NULL; // list is empty already
		E it = median(list, id, discrim); // Sort and take the median for "discrim" dimension
		// Erase "it"

		if (root == NULL)
		{
			root = new BSTNode<int, E>(0, it, NULL, NULL);
		}
		list.erase(list.begin() + id);
		size_t n = list.size();
		if (list.empty())
		{
			nodecount++;
			return root;
		}
		// Left list
		vector<E>leftList(list.begin(), list.begin() + n / 2);
		vector<E>rightList(list.begin() + n / 2, list.end());
		// Build left tree
		root->setLeft(balanceBuildHelp(root->left(), leftList, (discrim + 1) % D));
		// build right tree
		root->setRight(balanceBuildHelp(root->right(), rightList, (discrim + 1) % D));
		return root;
	}


	// Find the record with the given coordinates
	bool findhelp(BinNode<E>* root, int* coord,
		E& e, int discrim) const {
		// Member "coord" of a node is an integer array storing
		// the node's coordinates.
		if (root == NULL) return false;     // Empty tree
		int* currcoord = (root->element())->coord();
		if (EqualCoord(currcoord, coord)) { // Found it
			e = root->element();
			return true;
		}
		if (currcoord[discrim] < coord[discrim])
			return findhelp(root->left(), coord, e, (discrim + 1) % D);
		else
			return findhelp(root->right(), coord, e, (discrim + 1) % D);
	}

	// Query region w/int coordinates
	void queryRegionHelp(BinNode<E>* root, int* startCoord, int* endCoord, int discrim) const
	{
		if (root == NULL) return;
		int* currcoord = (root->element())->coord();
		// check if record at root inside the rectangle?
		if (InRect(currcoord, startCoord, endCoord))
			cout << root->element() << endl;
		if (currcoord[discrim] > endCoord[discrim])
			// don't need to check right child
			queryRegionHelp(root->left(), startCoord, endCoord, (discrim + 1) % D);
		else if (currcoord[discrim] < startCoord[discrim])
			queryRegionHelp(root->right(), startCoord, endCoord, (discrim + 1) % D);
		else
		{
			// Otherwise, we have to check both left and right children
			queryRegionHelp(root->left(), startCoord, endCoord, (discrim + 1) % D);
			queryRegionHelp(root->right(), startCoord, endCoord, (discrim + 1) % D);
		}
	}

	// Query region w/float coordinates
	void queryRegionHelp(BinNode<E>* root, float* startCoord, float* endCoord, int discrim) const
	{
		if (root == NULL) return;
		float* currcoord = (float*)((root->element())->coord());
		// check if record at root inside the rectangle?
		if (InRect(currcoord, startCoord, endCoord))
			cout << root->element() << endl;
		if (currcoord[discrim] > endCoord[discrim])
			// don't need to check right child
			queryRegionHelp(root->left(), startCoord, endCoord, (discrim + 1) % D);
		else if (currcoord[discrim] < startCoord[discrim])
			queryRegionHelp(root->right(), startCoord, endCoord, (discrim + 1) % D);
		else
		{
			// Otherwise, we have to check both left and right children
			queryRegionHelp(root->left(), startCoord, endCoord, (discrim + 1) % D);
			queryRegionHelp(root->right(), startCoord, endCoord, (discrim + 1) % D);
		}
	}

	// Query region w/result out as output and w/float coordinates
	void queryRegionGetResultHelp(BinNode<E>* root, float* startCoord, float* endCoord, int discrim, std::vector<E>& result) const
	{
		if (root == NULL) return;
		float* currcoord = (float*)((root->element())->coord());
		// check if record at root inside the rectangle?
		if (InRect(currcoord, startCoord, endCoord))
		{
			// If the record meets the thresholding criteria
			if (double(root->element()->oneMinusPercentile()) >= g_pFlatIdPt_thres)
				result.push_back(root->element()); // add to result when true

		}
		if (currcoord[discrim] > endCoord[discrim])
			// don't need to check right child
			queryRegionGetResultHelp(root->left(), startCoord, endCoord, (discrim + 1) % D, result);
		else if (currcoord[discrim] < startCoord[discrim])
			queryRegionGetResultHelp(root->right(), startCoord, endCoord, (discrim + 1) % D, result);
		else
		{
			// Otherwise, we have to check both left and right children
			queryRegionGetResultHelp(root->left(), startCoord, endCoord, (discrim + 1) % D, result);
			queryRegionGetResultHelp(root->right(), startCoord, endCoord, (discrim + 1) % D, result);
		}
	}

	// Print all points within distance "rad" of "coord"
	void regionhelp(BinNode<E>* root, int* coord,
		int rad, int discrim) const {
		if (root == NULL) return;       // Empty tree
		// Check if record at root is in circle
		if (InCircle((root->element())->coord(), coord, rad))
			cout << root->element() << endl;  // Do what is appropriate
		int* currcoord = (root->element())->coord();
		if (currcoord[discrim] > (coord[discrim] - rad))
			regionhelp(root->left(), coord, rad, (discrim + 1) % D);
		if (currcoord[discrim] < (coord[discrim] + rad))
			regionhelp(root->right(), coord, rad, (discrim + 1) % D);
	}
	// Return a pointer to the node with the least value in root
	// for the selected descriminator
	BinNode<E>* findmin(BinNode<E>* root,
		int discrim, int currdis) const {
		// discrim: discriminator key used for minimum search;
		// currdis: current level (mod D);
		if (root == NULL) return NULL;
		BinNode<E> *minnode = findmin(root->left(), discrim,
			(currdis + 1) % D);
		if (discrim != currdis) { // If not at descrim's level,
			// we must search both subtrees
			BinNode<E> *rightmin =
				findmin(root->right(), discrim, (currdis + 1) % D);
			// Check if right side has smaller key value
			minnode = min(minnode, rightmin, discrim);
		} // Now, minnode has the smallest value in children
		return min(minnode, root, discrim);
	}

	// Return in R the element (if any) which matches the coordinates K.
	// Return the updated subtree with R removed from the tree.
	BinNode<E>* removehelp(BinNode<E>* root, const int* K,
		BinNode<E>*& R, int descrim) {
		if (root == NULL) return NULL; // K is not in tree
		int* currcoord = (root->element())->coord();
		if (!EqualCoord(currcoord, K)) { // Continue searching
			if (K[descrim] < currcoord[descrim])
				root->setLeft(removehelp(root->left(), K, R, (descrim + 1) % D));
			else root->setRight(removehelp(root->right(), K, R, (descrim + 1) % D));
		}
		else {                            // Found K: remove this node
			R = root;
			if (root->right() == NULL) {    // No right child
				root->setRight(root->left()); //  so move left child to right
				root->setLeft(NULL);          // And make left NULL
			}
			if (root->right() == NULL)      // Both children are emtpy
				return NULL;                  // So the root will become NULL
			else {                          // Replace root with min on right
				BinNode<E>* temp;
				// First, find the min on right side
				temp = findmin(root->right(), descrim, (descrim + 1) % D);
				// Now, actually get that min node out of the right side
				root->setRight(removehelp(root->right(),
					(temp->element())->coord(), temp, (descrim + 1) % D));
				// At this point, we hold the minimum node, and the children are correct
				// So just replace root with the minimum and return the old root
				E te = root->element();
				root->setElement(temp->element());
				temp->setElement(te);
				R = temp;
			}
		}
		return root;
	}

	// Print out a KD tree in preorder
	void printhelp(BinNode<E>* root, int level) const {
		if (root == NULL) { // Empty tree
			cout << " NULL";
			return;
		}
		// Print node value
		cout << endl;
		for (int i = 0; i<level; i++)         // Indent to level
			cout << "  ";
		cout << root->element();

		// Print children
		printhelp(root->left(), level + 1);   // Do left subtree
		printhelp(root->right(), level + 1);  // Do right subtree
	}

	// Print to file
	void printToFileHelp(BinNode<E>* root, ofstream& ofFile, int level) const{
		if (root == NULL){ // empty tree
			ofFile << "NULL";
			return; 
		}
		// print node value
		ofFile << endl;
		for (int i = 0; i < level; i++)
			ofFile << " ";
		ofFile << root->element();
		// Print children
		printToFileHelp(root->left(), ofFile, level + 1);
		printToFileHelp(root->right(), ofFile, level + 1);
	}

};

#endif