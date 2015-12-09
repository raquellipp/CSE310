#include<list>
#include<iostream>
#include<stack>
#include<map>
#include<queue>
#include<vector>
#include<limits>
#include<iterator>
#include<stdio.h>

using namespace std;

//============================ Enumeration Mapping ============================
enum commandEnum{degree_distribution, shortest_path, diameter, components};
static std::map<std::string, commandEnum> mapCmdStringToEnum;
static void initializeCmdEnum(){
	mapCmdStringToEnum["degree-distribution"] = degree_distribution;
	mapCmdStringToEnum["shortest-path"] = shortest_path;
	mapCmdStringToEnum["diameter"] = diameter;
	mapCmdStringToEnum["components"] = components;
}
//=============================================================================

//============================== Global Variables =============================
// User Input -----------------------------------------------------------------
int vertexNum = 0; 				// Number of vertices
int edgeNum = 0;				// Number of edges

int commandCount = 0;			// Number of commands given by the user 
								// after the data has been inputted

string commandStr = "";			// The user inputted command		

// Linked-List ----------------------------------------------------------------
list<int> *adj; 				// List of pointers to linked-lists

// Degree ---------------------------------------------------------------------
int* vertexDegree;				// Stores the degree of each vertex

int* degreeCounts;				// Stores the number of times a degree (mapped 
								// by the index) occurs in vertexDegree

int maxDegree = -1;				// The highest degree in all of the nodes	

// Shortest-Path -------------------------------------------------------------
bool* discovered;				// Keeps track of the nodes that have been
								// discovered (used in shortest-path)

int *shortestDistance;			// Keeps track of the shortest distance (used
								// in shortest-path)

vector<int> includedVertex;		// Keeps track of all the vertices for which 
								// the shortest path has been computed.

vector<int> notIncludedVertex;	// Keeps track of all the vertices for which
								// the shortest path has not been computed
vector<int> shortestPathTrace;	// Keeps track of all the vertices for which
								// the shortest path has not been computed
int S = 0, T = 0;

// Components -----------------------------------------------------------------
enum color{white, gray, 
black};							// Color enumerations

struct vertexInfo{				// Defines a vertex struct that encompasses all
	int vertexData;				// information about each vertex into a single
	enum color vertexColor;		// object.
	struct vertexInfo *vertexPred;
	int vertexDistance;
};

queue<vertexInfo> vertexQueue;	// Vertex Queue

int connectionCount = 0; 
int connectionNum = 0;

// Diameter -------------------------------------------------------------------
int** distMatrix;
queue<int**> distMatrixQueue;
//=============================================================================

void BFS_intro(int, vector<vertexInfo>&);
void BFS(int, vector<vertexInfo>&);
void dijkstra(int, int);
void updateDistances();
void shortestPathTraceFunct(int);
void printAdj();

//=============================== Main Function ===============================
int main(){
//-----------------------------------------------------------------------------
// Store the number of vertices and the number of edges.
//-----------------------------------------------------------------------------	
	cin >> vertexNum;
	cin >> edgeNum;

//-----------------------------------------------------------------------------
// Define adj with a size of vertexNum. Then, use a for loop to store the
// information. Note that, because these graphs are undirected, each edge needs
// to be added to the linked list of both data values.
//-----------------------------------------------------------------------------
	adj = new list<int>[vertexNum]; // Remember to delete this later!

	// Store the information in adj
	int s, t;
	for (int index = 0; index < edgeNum; index++){
		cin >> s;
		cin >> t;
		adj[s].push_back(t);
		adj[t].push_back(s);
	}

//-----------------------------------------------------------------------------
// Define vertexDegree with a size of vertexNum. Then, use a for loop to store 
// the degree of each vertex in vertexDegree.
//-----------------------------------------------------------------------------	
	vertexDegree = new int[vertexNum]; // Remember to delete this later!
	for (int index = 0; index < vertexNum; index++){
		vertexDegree[index] = adj[index].size();
	}

//-----------------------------------------------------------------------------
// Initialize degreeCounts and set all indices to zero. This array will hold 
// the number of occurences of each degree that occures in vertexDegree. The
// index of the array will represent the degree, and the corresponding value
// will represent the number of occurences of that degree. 
//
// Note that a degree can be in the range [0, vertexNum]. So, degreeCounts 
// needs a size of vertexNum + 1. 
//-----------------------------------------------------------------------------
	degreeCounts = new int[vertexNum + 1]; // Remember to delete this later!
	for (int index = 0; index <= vertexNum; index++){
		degreeCounts[index] = 0; // Initializes the array to all zeros.
	}

	// Calculates the number of degrees and stores them in degreeCounts.
	for (int node = 0; node < vertexNum; node++){
		int degree = vertexDegree[node];
		degreeCounts[degree]++;
	}

//-----------------------------------------------------------------------------
// Find the max degree in degreeCounts. This will be used for the drawing of
// the histogram.
//-----------------------------------------------------------------------------	
	for (int degree = 0; degree <= vertexNum; degree++){
		if (degreeCounts[degree] > 0)
			maxDegree = degree;
	}


//-----------------------------------------------------------------------------
// Store the number of upcoming commands in the variable commandCount. Then,
// use a switch statement to implement these commands.
//-----------------------------------------------------------------------------	
	initializeCmdEnum();
	cin >> commandCount;
	cin >> commandStr;
	for (int cmd = 0; cmd < commandCount; cmd++){
		cout << endl << endl;
		switch(mapCmdStringToEnum[commandStr]){

//-----------------------------------------------------------------------------
// degree-distribution :
//
// Display the occurences of the degrees as a histogram. 
//
// The x-axis represents the degree, and the y-axis represents the number of 
// occurences of that degree.
//
// Draw the graph row-by-row (down the y-axis). Also note that the histogram
// only needs to be maxDegree units wide.
//-----------------------------------------------------------------------------
			case degree_distribution: {
				cout << "Degree (x) | Number of Occurrences (y)" << endl;
				for (int x = 1; x <= maxDegree; x++){
					if (degreeCounts[x] > 0){
						printf("%5d      | ", x);
						for (int i = 0; i < degreeCounts[x]; i++){
							cout << "*";
						}
						cout << endl;
					}
				}
				cout << endl;	
				break;
			}

//-----------------------------------------------------------------------------
// shortest-path s t :
//  
// s -> the source vertex
// t -> the destination vertex
// 
// Print the length of the shortest path in G between s and t, and also print
// the path to achieve it.
// 
// Dijkstra's single-source shortest path algorithm must be implemented to 
// compute the shortest path between s and t.
//
// v.d = shortest path between s and v.
// v.n = shortest path tree between s and v.
//-----------------------------------------------------------------------------	
			case shortest_path: {

				cout << "shortest-path was selected" << endl;
				
				// Input s and t.
				cin >> S;
				cin >> T;

				dijkstra(S,T);

				if (shortestDistance[T] < numeric_limits<int>::max() && 
					shortestDistance[T] > (-1)*(numeric_limits<int>::max())){
					cout << endl << "The shortest-path from " << S << " to " << T <<
						" is " << shortestDistance[T] << endl;
				
					shortestPathTrace.push_back(T);
					shortestPathTraceFunct(T);
					cout << endl << "Path:"; 
					for (int i=shortestPathTrace.size()-1; i>=0; i--){
						cout << " " << shortestPathTrace.at(i);
					}
				}

				else {
				 	cout << endl << "The shortest path from " << S << " to " << T <<
				 		" does not exist";
				}
					cout << endl;

				shortestPathTrace.clear();

				break;
			}

//-----------------------------------------------------------------------------
// components :
//  
// Uses the Breadth-first search to computed the number of connected components
// in the graph and their size.
//-----------------------------------------------------------------------------
			case components: {

				// Delete later!
				   vector<vertexInfo> vertexInfoList;
				for (int index = 0; index < vertexNum; index++){
					vertexInfo temp; // Delete later!
					vertexInfoList.push_back(temp); 
				}

				BFS_intro(0,vertexInfoList);
 			
				for (int index = 0; index < vertexNum; index++){
					if ((vertexInfoList[index].vertexColor) == white){
						BFS(index, vertexInfoList);
						if (connectionCount > 0){
							connectionNum++;
							cout << "Connection #" << connectionNum << " has "
								 << connectionCount << " connections." << endl;
						}
						connectionCount = 0;
					}
				}
				break;
			}

//-----------------------------------------------------------------------------
// diameter :
//
// Prints the diameter of the graph defined as the length of the longest 
// shortest path. The Floyd-Warshall all-pairs shortest path algorithm must be 
// used to compute the diameter; see §25.2 of our textbook for a discussion of 
// the Floyd-Warshall algorithm.
//-----------------------------------------------------------------------------
			case diameter: {
    			//All pair Shortest Path using Floyd Warshall's algorithm
    			distMatrix =  new int*[vertexNum];
    			for (int i=0; i<vertexNum; i++)
    				distMatrix[i] = new int[vertexNum];

   				 //initialize the distMatrix as Adjacency Matrix of the graph, i.e. 
    			//distMatrix[u][v] = 1 if (u,v) has a edge
    			//distMatrix[u][v] = 0 if u = v
    			//distMatrix[u][v] = numeric_limits<int>::max() if (u,v) does not have a edge

				for (int u=0; u<vertexNum; u++){
					for (int v=0; v<vertexNum; v++){
						distMatrix[u][v] = numeric_limits<int>::max();
					}
				}
    			for (int u=0; u<vertexNum; u++){
    				distMatrix[u][u] = 0;
    				for(list<int>::iterator v = adj[u].begin(); 
		 			v != adj[u].end(); ++v){
		 				distMatrix[u][*v] = 1;
		 				distMatrix[*v][u] = 1;
		 			}
    			}

    			distMatrixQueue.push(distMatrix);

    			//Algorithm Implementation
    			/* Add all vertices one by one to the set of intermediate vertices.
     			 ---> Before start of a iteration, we have shortest distances between all
      			pairs of vertices such that the shortest distances consider only the
      			vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      			----> After the end of a iteration, vertex no. k is added to the set of
     			 intermediate vertices and the set becomes {0, 1, 2, .. k} */
    			for (int k = 0; k < vertexNum; k++){
    				
    				int** curr = new int*[vertexNum];
    				curr = distMatrixQueue.back();
    				for (int i=0; i<vertexNum; i++)
    					curr[i] = distMatrixQueue.back()[i];
    				

    				int** prev = new int*[vertexNum];
    				prev = distMatrixQueue.back();
    				for (int i=0; i<vertexNum; i++)
    					prev[i] = distMatrixQueue.back()[i];


        			// Pick all vertices as source one by one
       				for (int i = 0; i < vertexNum; i++){
            			// Pick all vertices as destination for the
            			// above picked source
            			for (int j = 0; j < vertexNum; j++){
                			// If vertex k is on the shortest path from
                			// i to j, then update the value of dist[i][j] 

                			if (prev[i][k] == numeric_limits<int>::max() || 
                				prev[k][j] == numeric_limits<int>::max()){
                				curr[i][j] = prev[i][j];
                			}
                			else 
                				curr[i][j] = min(prev[i][j], (prev[i][k] + prev[k][j]));
            			}		
        			}

        			distMatrixQueue.push(curr);
   			 	}

   			 	// Find the longest shortest path
   			 	int** distMatrixN = new int*[vertexNum];
    			for (int i=0; i<vertexNum; i++)
    				distMatrixN[i] = new int[vertexNum];
    			distMatrixN = distMatrixQueue.back();

    			int diameterVal = distMatrixN[0][0];
    			for (int i=0; i<vertexNum; i++){
    				for (int j=0; j<vertexNum; j++){

    					if(distMatrixN[i][j] > diameterVal && (distMatrixN[i][j] != numeric_limits<int>::max()))
    						diameterVal = distMatrixN[i][j];
    				}
    			}

    			cout << endl << "The diameter is " << diameterVal << endl;
				break;
			}

			default: {
				cout << "Error. This is not a command." << endl;
				break;
			}
		}
		cin >> commandStr;
	}

	cout << endl << endl;

	delete[] adj;
	delete[] vertexDegree;
	delete[] degreeCounts;
	delete[] discovered;
	delete[] shortestDistance;
}
//=============================================================================


//============================== Other Functions ==============================
//-----------------------------------------------------------------------------
// BFS_intro :
// 
// Sets up vertexInfoList to get ready for BFS.
//-----------------------------------------------------------------------------
void BFS_intro(int S, vector<vertexInfo>& vertexInfoList2){
	connectionCount = 1;
	for (int vertex = 0; vertex < vertexNum; vertex++){
		if (vertex != S){
			vertexInfoList2.at(vertex).vertexData = vertex;
			vertexInfoList2.at(vertex).vertexColor = white;
			vertexInfoList2.at(vertex).vertexDistance = numeric_limits<int>::max();
			vertexInfoList2.at(vertex).vertexPred = NULL;
		}
	}
	vertexInfoList2.at(S).vertexData = S;
	vertexInfoList2.at(S).vertexColor = gray;
	vertexInfoList2.at(S).vertexDistance = numeric_limits<int>::max();
	vertexInfoList2.at(S).vertexPred = NULL;

	// Make sure the queue is clear
	while (!vertexQueue.empty()){
		vertexQueue.pop();
	}
}

//-----------------------------------------------------------------------------
// BFS :
// 
// Breadth-first search algorithm.
// 
// This algorithm works though the undirected graph to find connections
// between nodes. 
//
// Color meanings:
//	white -> Undiscovered
//	gray -> Has been searched once
//	black -> Has been searched twice
//-----------------------------------------------------------------------------
void BFS(int S, vector<vertexInfo>& vertexInfoList2){
	// Add S to the queue
	vertexQueue.push(vertexInfoList2.at(S));

	while (!vertexQueue.empty()){
		vertexInfo u = vertexQueue.front();
		vertexQueue.pop();

		 	for(list<int>::iterator it = adj[u.vertexData].begin(); 
		 		it != adj[u.vertexData].end(); ++it){
             	
             	if(vertexInfoList2.at(*it).vertexColor == white){
					vertexInfoList2.at(*it).vertexColor = gray;
					vertexInfoList2.at(*it).vertexDistance = 
						vertexInfoList2.at(u.vertexData).vertexDistance + 1;
					vertexInfoList2.at(*it).vertexPred = &u;
					vertexQueue.push(vertexInfoList2.at(*it));
				}
         	}

		vertexInfoList2.at(u.vertexData).vertexColor = black;
		connectionCount++;
	}
}

void dijkstra(int S, int T){
	// Initialized discovered to be all false.
	discovered = new bool[vertexNum]; // Delete this later!
	for(int index = 0; index < vertexNum; index++){
		discovered[index] = false;
	}
	// Set the start node to be discovered
	discovered[S] = true;

	// Initialize shortestDistance to be the max value of int
	// (infinity) for all values
	shortestDistance = new int[vertexNum]; // Delete this later!
	for (int index = 0; index < vertexNum; index++){
		shortestDistance[index] = numeric_limits<int>::max();
	}

	// Set the start node shortestDistance value to be 0
	// (since it is already discovered)
	shortestDistance[S] = 0;

	// Initialize includedVertex 
	includedVertex.clear();
	includedVertex.push_back(S);

	// Initialize notIncludedVertex
	notIncludedVertex.clear();
	for (int index = 0; index < vertexNum; index++){
		if (index != S)
			notIncludedVertex.push_back(index);
	}

	// *** DIJKSTRA'S SHORTEST PATH ALGORITHM ***
	while (notIncludedVertex.size() != 0){
		updateDistances();
		
		int u = notIncludedVertex.at(0);
		int indexSmall = 0;
		for (int i=1; i<notIncludedVertex.size(); i++){
			if (shortestDistance[u] > shortestDistance[notIncludedVertex.at(i)]){
				u = notIncludedVertex.at(i);
				indexSmall = i;
			}
		}

		discovered[u] = true;
		includedVertex.push_back(u);

		notIncludedVertex.erase(notIncludedVertex.begin() + indexSmall);

		if (u == T)
			break;
	}
}


void updateDistances(){
	for (int i = 0; i < includedVertex.size(); i++){
		for(list<int>::iterator it = adj[includedVertex.at(i)].begin(); 
		 		it != adj[includedVertex.at(i)].end(); ++it){
			if(!discovered[*it])
				if (shortestDistance[*it] > shortestDistance[includedVertex.at(i)] + 1)
					shortestDistance[*it] = shortestDistance[includedVertex.at(i)] + 1;
		}
	}
}	

void shortestPathTraceFunct(int current){

	while ((shortestDistance[current] != 0) && 
		(shortestDistance[current] <  numeric_limits<int>::max()) &&
		 (shortestDistance[current] > (-1)*numeric_limits<int>::max())){
		// new current starts here
		int temp = 0;
		for(list<int>::iterator it = adj[current].begin(); 
		it != adj[current].end(); ++it){
			if (shortestDistance[*it] == (shortestDistance[current] - 1)){
				shortestPathTrace.push_back(*it);
				temp = *it;
			}
		}
		current = temp;
	}
}

void printAdj(){
	cout << endl << endl << "adj contents:" << endl;
	for (int i = 0; i < vertexNum; i++){
		cout << i << ":";
		for(list<int>::iterator it = adj[i].begin(); 
		it != adj[i].end(); ++it){
			cout << " " << *it;
		}
		cout << endl;
	}
}

//=============================================================================
