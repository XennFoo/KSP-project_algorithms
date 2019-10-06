#include "pch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <algorithm>    
#include <vector>
#include <stdlib.h>
#include <math.h>
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
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/copy.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include "seq/knheap.c"
#include "bin/bheap.h"
#pragma warning(disable : 4996)
#define HTYPE KNHeap<double, vertex_descriptor> 
#define HINIT heap(INT_MAX, -INT_MAX)


typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::no_property, EdgeWeightProperty > DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<DirectedGraph>::edge_iterator edge_iterator;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<DirectedGraph>::out_edge_iterator out_edge_iterator;
DirectedGraph g;
DirectedGraph reverse_g;
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

std::vector<vertex_descriptor> generatePathFromAncestors2(vertex_descriptor s, vertex_descriptor f, std::vector<vertex_descriptor> p) {
	std::vector<vertex_descriptor> path;
	path.push_back(f);
	while (f != s) {
		f = p[f];
		path.push_back(f);
	}
	return path;
}


vertex_descriptor yellowDijkstra(DirectedGraph& g, vertex_descriptor s, vertex_descriptor f, std::vector<std::vector<vertex_descriptor>>& express, std::vector<vertex_descriptor>& parent, std::vector<double>& distance, std::vector<double>& reald, std::vector<int>& color, std::vector<double> reverse_d, HTYPE& heap) {
	vertex_iterator vi, vend;
	heap.clear();
	int flag = 0;
	vertex_descriptor temp2, temp, result;
	temp = NULL;
	boost::graph_traits<DirectedGraph>::out_edge_iterator ei, ei_end;
	//initializing distance map

	for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) {
		distance[*vi] = INT_MAX;
		reald[*vi] = INT_MAX;
		parent[*vi] = NULL;

	}
	distance[s] = 0;
	reald[s] = 0;
	heap.insert(0, s);
	while (heap.getSize() != 0) {
		double k;
		vertex_descriptor v;
		heap.deleteMin(&k, &v);
		temp = v;
		//check if we reached a final edge
		for (auto& row : express) {
			for (auto& col : row) {
				if (temp == col && temp != s) {
					result = col;
					flag = 1;
				}
				if (temp == f) {
					result = f;
					flag = 1;
				}
			}
		}
		if (flag) {
			break;
		}
		boost::tie(ei, ei_end) = boost::out_edges(temp, g);
		for (; ei != ei_end; ++ei) {
			temp2 = boost::target(*ei, g);
			edge_descriptor e = edge(temp, temp2, g).first;
			//the weight is considered originalweight + d(head) - d(tail)
			double weight = get(boost::edge_weight_t(), g, e) + reverse_d[temp2] - reverse_d[temp];
			double realweight = get(boost::edge_weight_t(), g, e);
			if (distance[temp] + weight < distance[temp2] && color[temp2] == 1) {
				distance[temp2] = distance[temp] + weight;
				reald[temp2] = reald[temp] + realweight;
				parent[temp2] = temp;
				heap.insert(distance[temp2], temp2);
			}
		}
	}
	return result;
};


int main(int, char *[])
{
	//read graph from DIMACS file
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


	//FENGS ALGORITHM
	boost::copy_graph(boost::make_reverse_graph(g), reverse_g);
	std::cout << "If you want to load destinations from a file press (2), otherwise press (1) " << std::endl;
	int choice, Numberoftimes, S, F, K;
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
		std::vector<vertex_descriptor> path;
		std::vector< std::vector<vertex_descriptor>> ListA; //the list that contains the k best paths
		std::vector<routes> paths; // this contains the k best paths except the first one as routes 
		std::priority_queue<routes, std::vector<routes>, comperator> ListB; //the list that contains all the candidate paths
		boost::graph_traits < DirectedGraph >::vertex_iterator vi, vend;



		//finds the first path
		dijkstra_shortest_paths(g, s,
			predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
			distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
		path = generatePathFromAncestors(s, f, p);
		ListA.push_back(path);


		//reversing the graph in order to find T
		std::vector<vertex_descriptor> reverse_p(num_vertices(g));
		std::vector<double> reverse_d(num_vertices(g));
		dijkstra_shortest_paths(reverse_g, f,
			predecessor_map(boost::make_iterator_property_map(reverse_p.begin(), get(boost::vertex_index, g))).
			distance_map(boost::make_iterator_property_map(reverse_d.begin(), get(boost::vertex_index, g))));
		boost::graph_traits < DirectedGraph >::edge_iterator ei, eend;
		vertex_descriptor v;

		//finding the childrens tree
		std::vector<std::vector<vertex_descriptor>> children_tree(num_vertices(g));
		boost::make_iterator_property_map(children_tree.begin(), get(boost::vertex_index, g));
		vertex_iterator vs, vsend;
		for (boost::tie(vs, vsend) = vertices(g); vs != vsend; ++vs) {
			if (reverse_p[*vs] != NULL) {
				v = reverse_p[*vs];
				children_tree[v].push_back(*vs);//now the children of every node are easily accesible
			}

		}


		//iteration to find the k paths(find candidate paths)
		for (int i = 0; i < K - 1; i++) {
			//findCandidates from path Pi

			//find the distance from s to every node of the path

			path.clear();
			path = ListA[i];
			std::vector<double> prevDistance(path.size());
			std::vector<double> prevDistanceV(num_vertices(g));
			boost::make_iterator_property_map(prevDistanceV.begin(), get(boost::vertex_index, g));
			prevDistance[0] = 0;
			for (int j = 1; j < path.size(); j++) {
				std::pair<edge_descriptor, bool> ed = boost::edge(path[j - 1], path[j], g);
				double weight = get(boost::edge_weight_t(), g, ed.first);
				prevDistance[j] = prevDistance[j - 1] + weight;
			}

			int j = 0;
			for (auto it = begin(path); it != end(path); ++it) {
				prevDistanceV[*it] = prevDistance[j];
				j++;
			}

			//create a color_map and paint every node green
			std::vector<int> color(num_vertices(g));
			boost::make_iterator_property_map(color.begin(), get(boost::vertex_index, g));
			std::vector<std::vector<vertex_descriptor>> express(num_vertices(g));
			boost::make_iterator_property_map(express.begin(), get(boost::vertex_index, g));
			vertex_iterator v, vend;
			for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
				color[*v] = 0; // 0 is the green color  yellow is 1   red is 2
			}
			std::vector<vertex_descriptor> tempath;
			vertex_descriptor u;
			if (i != 0) { //if this isn't the first path we need to remove edges
				u = paths[i - 1].getDevNode();


				//  ###  findEdgesToRemove  ###
				//removes all the edges that give us a path that allready exists
				std::vector<vertex_descriptor> rem = paths[i - 1].getRemPath();
				std::pair<edge_descriptor, bool> ed = boost::edge(rem[0], rem[1], g);
				double weight2 = get(boost::edge_weight_t(), g, ed.first);
				boost::put(boost::edge_weight_t(), g, ed.first, weight2 + INT_MAX / 2);
				int k = paths[i - 1].getKpath();//******************************/
				if (k > 0) {
					routes father = paths[k - 1];
					while (u == father.getDevNode()) {
						rem = father.getRemPath();
						//cutting the edge
						ed = boost::edge(rem[0], rem[1], g);
						weight2 = get(boost::edge_weight_t(), g, ed.first);
						boost::put(boost::edge_weight_t(), g, ed.first, weight2 + INT_MAX / 2);
						k = father.getKpath();
						if (k == 0) {
							break;
						}
						father = paths[k - 1];
					}
				}
				int countc = 0;
				std::vector<vertex_descriptor> cpath = ListA[k];
				for (auto it = begin(cpath); *it != u; ++it) {
					countc++;
				}
				ed = boost::edge(cpath[countc], cpath[countc + 1], g);
				weight2 = get(boost::edge_weight_t(), g, ed.first);
				boost::put(boost::edge_weight_t(), g, ed.first, weight2 + INT_MAX / 2);
			}
			else {
				u = s;
			}

			std::vector<vertex_descriptor> redpath;
			if (i != 0) {//if it its not the first path, paint all the nodes until node u red
				int kpath = paths[i - 1].getKpath();
				redpath = ListA[kpath];
			}
			else {
				redpath.push_back(s);
			}
			std::vector<vertex_descriptor> S;
			std::vector<vertex_descriptor> Y;

			for (auto it = begin(redpath); *it != u; ++it) {
				color[*it] = 2;  //change the paths nodes to red
			}
			color[u] = 2; // paint the vertex u red


			for (auto it = begin(redpath); *it != u; ++it) {
				//color[*it] = 2;  //change the paths nodes to red



				//  ###  getUpstreamNodes();  ###  
				//here I paint every green upstream node yellow
				S.push_back(*it);
				while (!S.empty()) {
					vertex_descriptor t;
					t = S.back();
					S.pop_back();
					for (auto itt = begin(children_tree[t]); itt != end(children_tree[t]); ++itt) {
						if (color[*itt] == 0) {
							S.push_back(*itt);
							Y.push_back(*itt);
						}
					}
				}
			}
			S.push_back(u);
			while (!S.empty()) {
				vertex_descriptor t;
				t = S.back();
				S.pop_back();
				for (auto itt = begin(children_tree[t]); itt != end(children_tree[t]); ++itt) {
					if (color[*itt] == 0) {
						S.push_back(*itt);
						Y.push_back(*itt);
					}
				}
			}
			for (auto it = begin(Y); it != end(Y); ++it) {
				color[*it] = 1;

			}


			//  ###  findExpressEdges();  ###
			//transform express edges to final a final edge
			out_edge_iterator eii, eii_end;
			vertex_descriptor t;
			for (auto it = begin(Y); it != end(Y); ++it) {
				for (boost::tie(eii, eii_end) = out_edges(*it, g); eii != eii_end; ++eii) {
					t = boost::target(*eii, g);
					if (color[t] == 0) {
						express[*it].push_back(t);
						color[t] = 1;
					}
				}
			}
			for (boost::tie(eii, eii_end) = out_edges(u, g); eii != eii_end; ++eii) {
				t = boost::target(*eii, g);
				if (color[t] == 0) {
					express[u].push_back(t);
					color[t] = 1;
				}
			}
			color[f] = 1;
			std::vector<vertex_descriptor> pathToCheck;
			if (i == 0) {
				pathToCheck = ListA[0];
			}
			else {
				pathToCheck = paths[i - 1].getRemPath();
			}
			int vertex_counter = 0;
			int whilecounter = 0;
			while (true) {
				//remove the edge u->v
				std::pair<edge_descriptor, bool> ed2 = boost::edge(pathToCheck[vertex_counter], pathToCheck[vertex_counter + 1], g);
				double weight = get(boost::edge_weight_t(), g, ed2.first);
				if (weight < INT_MAX / 2) {
					boost::put(boost::edge_weight_t(), g, ed2.first, weight + INT_MAX / 2);
				}

				//call Dijkstra from point u(pathToCheck[vertex_counter]) to some endPoint (express)
				//get the resulting path into ListB(routes)
				std::vector<vertex_descriptor> pp(num_vertices(g));
				std::vector<double> dd(num_vertices(g));
				std::vector<double> reald(num_vertices(g));
				boost::make_iterator_property_map(pp.begin(), get(boost::vertex_index, g));
				boost::make_iterator_property_map(dd.begin(), get(boost::vertex_index, g));
				boost::make_iterator_property_map(reald.begin(), get(boost::vertex_index, g));
				vertex_descriptor finalpoint;
				for (int g = 0; g < vertex_counter + 1; g++) {
					color[pathToCheck[g]] = 2;

				}
				finalpoint = yellowDijkstra(g, pathToCheck[vertex_counter], f, express, pp, dd, reald, color, reverse_d, heap);

				std::vector<vertex_descriptor> path2 = generatePathFromAncestors(pathToCheck[vertex_counter], finalpoint, pp);
				routes newroute = routes(i, pathToCheck[vertex_counter], path2, reald[finalpoint] + prevDistanceV[pathToCheck[vertex_counter]] + reverse_d[finalpoint]);
				ListB.push(newroute);
				pp.clear();
				dd.clear();
				heap.clear();
				if (pathToCheck[vertex_counter + 1] == f) {
					break;

				}
				u = pathToCheck[vertex_counter + 1];
				color[u] = 2;

				//  ###  GetUpstreamNodes  ### for vertex u 
				Y.clear();
				S.push_back(u);
				while (!S.empty()) {
					vertex_descriptor t;
					t = S.back();
					S.pop_back();
					for (auto itt = begin(children_tree[t]); itt != end(children_tree[t]); ++itt) {
						if (color[*itt] == 0) {
							S.push_back(*itt);
							Y.push_back(*itt);
						}
					}
				}
				for (auto it = begin(Y); it != end(Y); ++it) {
					color[*it] = 1;

				}
				//recover all express edges in eXpress v e Y
				for (auto it = begin(Y); it != end(Y); ++it) {
					if (!express[*it].empty()) {
						//paint back the express nodes green
						for (auto itt = begin(express[*it]); itt != end(express[*it]); ++itt) {
							if (color[*itt] == 1) {
								color[*itt] = 0;
							}
						}
						express[*it].clear();
					}
				}

				//  ###  FindExpressEdges ### for Yu + u
				for (auto it = begin(Y); it != end(Y); ++it) {

					for (boost::tie(eii, eii_end) = out_edges(*it, g); eii != eii_end; ++eii) {
						t = boost::target(*eii, g);
						if (color[t] == 0) {
							express[*it].push_back(t);
							color[t] = 1;
						}

					}
				}
				for (boost::tie(eii, eii_end) = out_edges(u, g); eii != eii_end; ++eii) {
					t = boost::target(*eii, g);
					if (color[t] == 0) {
						express[u].push_back(t);
						color[t] = 1;
					}
				}
				vertex_counter++;
				whilecounter++;
			}
			//Getting the k path in ListA
			routes temporary = ListB.top();
			ListB.pop();
			for (boost::tie(ei, eend) = boost::edges(g); ei != eend; ++ei)
			{
				double weight = get(boost::edge_weight_t(), g, *ei);
				if (weight > INT_MAX / 2) {
					boost::put(boost::edge_weight_t(), g, *ei, weight - INT_MAX / 2);
				}
			}

			color.clear();
			express.clear();


			std::vector<vertex_descriptor> restpath = temporary.getRemPath();
			vertex_descriptor finalpoint = restpath.back();
			//we get the rest of the path from the reversed graph
			std::vector<vertex_descriptor> remainingPath = generatePathFromAncestors2(f, finalpoint, reverse_p);
			int kpath = temporary.getKpath();
			vertex_descriptor dnode = temporary.getDevNode();
			std::vector<vertex_descriptor> prevpath;
			prevpath = ListA[kpath];
			std::vector<vertex_descriptor> newPath, newPathonPaths;

			//find the new path by connecting Qa and Qb
			for (auto it = begin(prevpath); *it != dnode; ++it) {
				newPath.push_back(*it);
			}
			for (auto z = begin(restpath); z != end(restpath); ++z) {
				newPath.push_back(*z);
				newPathonPaths.push_back(*z);

			}
			newPath.pop_back(); //this erases the finalPoint that is both the end of restpath and the start of remainingPath
			newPathonPaths.pop_back();
			for (auto it = begin(remainingPath); it != end(remainingPath); ++it) {
				newPath.push_back(*it);
				newPathonPaths.push_back(*it);
			}

			routes entry = routes(temporary.getKpath(), temporary.getDevNode(), newPathonPaths, temporary.getDistance());
			paths.push_back(entry);
			ListA.push_back(newPath);
		}
		//Prints the path that gets confirmed as the i+1 shortest path
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
		}*/
	}

	auto end = std::chrono::high_resolution_clock::now();
	// Calculating total time taken by the program. 
	double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	std::cout << "choice is :" << choice << std::endl;
	std::cout << "average time taken by Fengs's algorithm : " << std::fixed
		<< time_taken / choice << std::setprecision(9);
	std::cout << " sec" << std::endl;
	std::cout << "Programm exits";
	return 0;
}
