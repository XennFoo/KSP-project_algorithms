
#include "pch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <algorithm>    
#include <vector>
#include <chrono> 
#include<boost/tokenizer.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include "seq/knheap.c"
#include "bin/bheap.h"
#define HTYPE KNHeap <double, vertex_descriptor> 
#define HINIT heap(INT_MAX, -INT_MAX)
#pragma warning(disable : 4996)


typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::no_property, EdgeWeightProperty > DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<DirectedGraph>::edge_iterator edge_iterator;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor edge_descriptor;
DirectedGraph g;

std::vector<vertex_descriptor> generatePathFromAncestors(vertex_descriptor s, vertex_descriptor f, std::vector<vertex_descriptor> p);

//ListB objects
class routes {
	int kpath;
	vertex_descriptor devnode;
	std::vector<vertex_descriptor> rempath;
	double distance;
public:
	routes(int _kpath, vertex_descriptor _devnode, std::vector<vertex_descriptor> _rempath, double _distance) {
		kpath = _kpath;
		devnode = _devnode;
		rempath = _rempath;
		distance = _distance;
	}
	int getKpath() const { return kpath; }
	vertex_descriptor getDevNode() const { return devnode; }
	std::vector<vertex_descriptor> getRemPath() const { return rempath; }
	double getDistance() const { return distance; }
};

//Helps the operation of minheap
class comperator {
public:
	int operator()(const routes& r1, const routes& r2) {
		return r1.getDistance() > r2.getDistance();
	}
};

//returns a path from a parent tree
std::vector<vertex_descriptor> generatePathFromAncestors(vertex_descriptor s, vertex_descriptor f, std::vector<vertex_descriptor> p) {
	std::vector<vertex_descriptor> reverse_path;
	reverse_path.push_back(f);
	while (f != s) {
		f = p[f];
		reverse_path.push_back(f);
	}
	std::reverse(reverse_path.begin(), reverse_path.end());
	return reverse_path;
}


void myDijkstra(DirectedGraph& g, vertex_descriptor s, vertex_descriptor f, std::vector<vertex_descriptor>& parent, std::vector<double>& distance, std::vector<int>& cyclecheck, HTYPE& heap) {
	vertex_iterator vi, vend;

	heap.clear();
	vertex_descriptor temp2, temp;
	temp = NULL;
	boost::graph_traits<DirectedGraph>::out_edge_iterator ei, ei_end;
	//initializing distance map

	for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
		distance[*vi] = INT_MAX;
		parent[*vi] = NULL;

	}
	distance[s] = 0;
	heap.insert(0, s);
	while (heap.getSize() != 0) {
		double k;
		vertex_descriptor v;
		heap.deleteMin(&k, &v);
		temp = v;
		if (temp == f) {
			break;
		}
		boost::tie(ei, ei_end) = boost::out_edges(temp, g);
		for (; ei != ei_end; ++ei) {
			temp2 = boost::target(*ei, g);
			edge_descriptor e = edge(temp, temp2, g).first;
			double weight = get(boost::edge_weight_t(), g, e);
			if (distance[temp] + weight < distance[temp2] && cyclecheck[temp2]) {
				distance[temp2] = distance[temp] + weight;
				parent[temp2] = temp;
				heap.insert(distance[temp2], temp2);
			}
		}
	}
	return;
};


int main(int, char *[])
{
	//readgraph from DIMACS file
	std::ifstream myfile;
	char filename[50];
	std::cout << "Enter the file name: ";
	std::cin.getline(filename, 50);
	myfile.open(filename);
	std::string line;
	std::cout << "Loading graph structure....";
	while (std::getline(myfile, line)) {
		std::istringstream iss(line);
		boost::tokenizer<> tok(line);
		boost::tokenizer<>::iterator beg = tok.begin();
		if (*beg == "a") {
			++beg;
			int from = std::stoi(*beg);
			++beg;
			int to = std::stoi(*beg);
			++beg;
			int weight = std::stoi(*beg);
			boost::add_edge(from, to, weight, g);
		}
	}
	std::cout << "Loading finished!" << std::endl;
	HTYPE HINIT;

	//YENS ALGORITHM

	std::cout << "If you want to load destinations from a file press (2), otherwise press (1) " << std::endl;
	int choice, S, F, K;
	scanf("%d", &choice);
	while (choice != 1 && choice != 2) {
		std::cout << " wrong choice, try again" << std::endl;
		scanf("%d", &choice);
	}
	std::ifstream myfile2;
	if (choice == 2) {
		myfile2.open("testfile.txt");
		std::getline(myfile2, line);
		std::istringstream iss(line);
		boost::tokenizer<> tok(line);
		boost::tokenizer<>::iterator beg = tok.begin();
		if (*beg == "n") {
			++beg;
			choice = std::stoi(*beg);
		}
		std::getline(myfile2, line);
		std::istringstream iss2(line);
		boost::tokenizer<> tok2(line);
		beg = tok2.begin();
		if (*beg == "k") {
			++beg;
			K = std::stoi(*beg);
		}
	}
	else {
		std::cout << "give the number k of the minimum paths   ";
		scanf("%d", &K);
	}
	auto start = std::chrono::high_resolution_clock::now();
	for (int testtimes = 0; testtimes < choice; testtimes++) {
		if (choice == 1) {
			std::cout << "Give the start node   ";
			scanf("%d", &S);
			std::cout << "Give the final destination   ";
			scanf("%d", &F);
		}
		else {
			std::getline(myfile2, line);
			std::istringstream iss(line);
			boost::tokenizer<> tok(line);
			boost::tokenizer<>::iterator beg = tok.begin();
			if (*beg == "p") {
				++beg;
				S = std::stoi(*beg);
				++beg;
				F = std::stoi(*beg);
				std::cout << "searching from " << S << " to " << F << std::endl;
			}
		}
		std::vector<vertex_descriptor> p(num_vertices(g));
		std::vector<double> d(num_vertices(g));
		vertex_descriptor s = vertex(S, g);
		vertex_descriptor f = vertex(F, g);
		dijkstra_shortest_paths(g, s,
			predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
			distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
		std::vector<vertex_descriptor> path = generatePathFromAncestors(s, f, p);
		std::vector< std::vector<vertex_descriptor>> ListA;
		std::vector<vertex_descriptor> current_path;
		std::vector<routes> ListofRoutes;
		std::priority_queue<routes, std::vector<routes>, comperator> ListB;
		ListA.push_back(path);
		routes temproute = routes(-1, s, path, d[f]);
		ListofRoutes.push_back(temproute);
		boost::graph_traits < DirectedGraph >::vertex_iterator vi, vend;

		//Iteration for every path to be found
		for (int k = 2; k <= K; k++) {
			path.clear();
			path = ListA.back();  //k-1 path
			routes routepath = ListofRoutes.back();
			int startingpoint = 0;
			std::vector<int> cyclecheck(num_vertices(g));
			boost::make_iterator_property_map(cyclecheck.begin(), get(boost::vertex_index, g));
			vertex_iterator v, vend;
			for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
				cyclecheck[*v] = 1; 
			}
			//stores each nodes distance from the root
			std::vector<double> prevDistance(path.size());
			prevDistance[0] = 0;
			for (int i = 1; i < (int)path.size();i++) {
				prevDistance[i] = 0;
				std::pair<edge_descriptor, bool> ed = boost::edge(path[i - 1], path[i], g);
				double weight = get(boost::edge_weight_t(), g, ed.first);
				prevDistance[i] = prevDistance[i - 1] + weight;
				if (path[i - 1] == routepath.getDevNode()) {
					startingpoint = i - 1;
				}

			}
			vertex_descriptor deviat_node = routepath.getDevNode();
			if (k > 2) { //if this isn't the first path we need to remove edges
				vertex_descriptor u = routepath.getDevNode();

				//  ###  findEdgesToRemove  ###
				//removes all the edges that give us a path that allready exists
				std::vector<vertex_descriptor> rem = routepath.getRemPath();
				std::pair<edge_descriptor, bool> ed = boost::edge(rem[0], rem[1], g);
				double weight2 = get(boost::edge_weight_t(), g, ed.first);
				boost::put(boost::edge_weight_t(), g, ed.first, weight2 + INT_MAX / 2);
				int kk = routepath.getKpath();//******************************/
				if (kk > -1) {
					routes father = ListofRoutes[kk];
					while (u == father.getDevNode() && k != -1) {
						rem = father.getRemPath();
						//cutting the edge
						ed = boost::edge(rem[0], rem[1], g);
						if (weight2 < INT_MAX / 2) {
							boost::put(boost::edge_weight_t(), g, ed.first, weight2 + INT_MAX / 2);
						}
						kk = father.getKpath();
						if (kk == -1) {
							break;
						}
						father = ListofRoutes[kk];
					}
					int countc = 0;
					if (kk != -1) {
						std::vector<vertex_descriptor> cpath = ListA[kk];
						for (auto it = begin(cpath); *it != u; ++it) {
							countc++;
						}
						ed = boost::edge(cpath[countc], cpath[countc + 1], g);
						weight2 = get(boost::edge_weight_t(), g, ed.first);
						if (weight2 < INT_MAX / 2) {
							boost::put(boost::edge_weight_t(), g, ed.first, weight2 + INT_MAX / 2);
						}
					}
				}
				//marks the nodes of the path that shouldnt be accessed
				for (int j = 0; j < startingpoint; j++) {
					cyclecheck[path[j]] = 0;
				}
			}
			//for every node of the last path (exept the last one)
			//nodes before the deviation node have already been examined so they are skipped
			for (int i = startingpoint; i < (int)path.size() - 1; i++) {
				std::pair<edge_descriptor, bool> ed = boost::edge(path[i], path[i + 1], g);
				double weight = get(boost::edge_weight_t(), g, ed.first);
				if (weight < INT_MAX / 2) {
					boost::put(boost::edge_weight_t(), g, ed.first, weight + INT_MAX / 2);
				}


				//PART B OF ALGORITHM
				//for every node of the previous path find the new best route and store it to ListB
				vertex_descriptor sstart = path[i];
				cyclecheck[sstart] = 0;
				std::vector<vertex_descriptor> p2(num_vertices(g));
				std::vector<double> d2(num_vertices(g));
				boost::make_iterator_property_map(p2.begin(), get(boost::vertex_index, g));
				boost::make_iterator_property_map(d2.begin(), get(boost::vertex_index, g));
				myDijkstra(g, sstart, f, p2, d2, cyclecheck, heap);
				std::vector<vertex_descriptor> path2 = generatePathFromAncestors(sstart, f, p2);
				routes newroute = routes(k - 2, sstart, path2, d2[f] + prevDistance[i]);
				ListB.push(newroute);

			}


			//Reverting edge weights
			boost::graph_traits < DirectedGraph >::edge_iterator ei, eend;
			for (boost::tie(ei, eend) = boost::edges(g); ei != eend; ++ei)
			{
				double weight = get(boost::edge_weight_t(), g, *ei);
				if (weight > INT_MAX / 2) {
					boost::put(boost::edge_weight_t(), g, *ei, weight - INT_MAX / 2);
				}
			}

			//Getting the k path in ListA
			routes temporary = ListB.top();
			ListB.pop();
			std::vector<vertex_descriptor> restpath = temporary.getRemPath();
			int kpath = temporary.getKpath();
			vertex_descriptor dnode = temporary.getDevNode();
			path.clear();
			path = ListA[kpath];
			std::vector<vertex_descriptor> newPath;

			//find the new path by connecting Qa and Qb
			for (auto it = begin(path); *it != dnode; ++it) {
				newPath.push_back(*it);
			}
			for (auto z = begin(restpath); z != end(restpath); ++z) {
				newPath.push_back(*z);

			}
			ListA.push_back(newPath); //put the new path in ListA
			ListofRoutes.push_back(temporary);
			//std::cout << "Finished for k = " << k << std::endl;
		}
		/*int z = 0;
		while (z < K) {
			std::cout << z + 1 << " path" << std::endl;
			path.clear();
			path = ListA[z];
			for (auto it = begin(path); it != end(path); ++it) {
				std::cout << "-->" << *it;
			}
			z++;
			std::cout << std::endl;
		}
		*/
	}
	auto end = std::chrono::high_resolution_clock::now();
	// Calculating total time taken by the program. 
	double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	std::cout << "average time taken by Yen's algorithm : " << std::fixed
		<< time_taken/choice << std::setprecision(9);
	std::cout << " sec" << std::endl;
	std::cout << "Programm exits";
	return 0;
}