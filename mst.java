package comp3600_ass_2;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class mst {
	public static void main(String[] args) {
		runMSTTest(10);
		runMSTTest(100);
		runMSTTest(500);
		runMSTTest(1000);
	}
	
	private static void runMSTTest(int n) {
		double krusAveTime = 0;
		double krusAveL = 0;
		double primAveTime = 0;
		double primAveL= 0;
		double[][][] graphs = new double[20][n][n];
		
		//generate all 20 graphs, stored into an array
		for (int i = 0; i < 20; i++) {
			graphs[i] = initGraph(n);
		}

		for (int i = 0; i < 20; i++) {
			long startTime = 0;
			long endTime = 0;
			
			//run Kruskal's and Prim's algorithms on the graph, returning average weight L(n) and time to complete
			startTime = System.currentTimeMillis();
			krusAveL += kruskal(graphs[i]);
			
			/*if(n == 1000) {
				System.out.println("KrusAveSum: "+krusAveL);
			}*/
			
			endTime = System.currentTimeMillis();
			
			krusAveTime += (endTime - startTime);
			
			startTime = System.currentTimeMillis();
			//primAveL += prim(graphs[i]);
			
			/*if(n == 1000) {
				System.out.println("PrimAveSum: "+primAveL);
			}*/
			
			endTime = System.currentTimeMillis();
			
			primAveTime += (endTime - startTime);
//			primL += primVals[1];
		}
		
		//Calculate averages
		krusAveTime = krusAveTime / 20;
		krusAveL = krusAveL / 20;
		//primAveTime = primAveTime / 20;
		//primAveL = primAveL / 20;
		
		System.out.println("n = " + String.valueOf(n) + " --> Kruskal Ave Time: " + krusAveTime + ", Kruskal Ave L(n): " + krusAveL + "   |||   Prim Ave Time: " + primAveTime + ", Prim Ave L(n): " + primAveL);
	}
	
	//generates a complete symmetric graph, randomly weighted between 0 and 1
	private static double[][] initGraph(int n) {
		double[][] graph = new double[n][n];
		
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (i != j) {
					graph[i][j] = Math.random();
					graph[j][i] = graph[i][j];
				}
			}
		}
		return graph;
	}
	
	//prints out a given graph matrix
	private static void printG(double[][] g, int n) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				System.out.print(g[i][j] + "  ");
				
			}
			System.out.print("\n");
		}
	}
	
	
	//KRUSKAL'S ALGORITHM
	private static double kruskal(double[][] g) {
		int n = g[0].length;
		double[][] mst = new double[n][n];
		DisjointSet[] forest = new DisjointSet[n];
		int numEdges = (int) ((Math.pow(n,2)) - n)/2;
		Edge[] edges = new Edge[numEdges];
		double krusAveL = 0;
		
		//create DisjointSet for each vertex
		for (int i = 0; i < n; i++) {
			DisjointSet t = new DisjointSet(i);
			forest[i] = t;
		}
		
		//create array of edges
		int nextEdge = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				edges[nextEdge] = new Edge(g[i][j], i, j);
				nextEdge++;
			}
		}
		
		//sort edges
		edges = quickSort(edges, 0, numEdges-1);
		
		//run Kruskal's
		for(int i = 0; i < numEdges; i++) {
			DisjointSet treeSrc = null;
			DisjointSet treeDest = null;
			
			/*for(DisjointSet t : forest) {
				if(t.value == edges[i].src) {
					treeSrc = findSet(t);
				}
				if(t.value == edges[i].dest) {
					treeDest = findSet(t);
				}
				
				if((treeSrc != null) && (treeDest != null)) {
					if(treeSrc != treeDest) {
						//add edge to MST
						mst[edges[i].src][edges[i].dest] = edges[i].weight;
						mst[edges[i].dest][edges[i].src] = edges[i].weight;
						
						krusAveL += edges[i].weight;
						
						//union sets
						union(treeSrc, treeDest);
						//forest.remove(treeDest);
						break;
					}else break;
				}
			}*/
			treeSrc = findSet(forest[edges[i].src]);
			treeDest = findSet(forest[edges[i].dest]);
			
			if((treeSrc != null) && (treeDest != null)) {
				if(treeSrc != treeDest) {
					//add edge to MST
					mst[edges[i].src][edges[i].dest] = edges[i].weight;
					mst[edges[i].dest][edges[i].src] = edges[i].weight;
					
					krusAveL += edges[i].weight;
					
					/*if(n == 1000) {
						System.out.println("Edge weight "+ i + ": " + edges[i].weight);
					}*/
					
					//union sets
					union(treeSrc, treeDest);
					//forest.remove(treeDest);
				}
			}
		}
		//long endTime = System.currentTimeMillis();

		//krusAveL = krusAveL / (n-1); //n-1 == number of edges in MST

		return krusAveL;
	}
	
	
//	//PRIM'S ALGORITHM
//	private static double prim(double[][] g) {
//		int n = g[0].length;
//		double[][] mst = new double[n][n];
//		//DisjointSet tree = null;
//		int[] tree = new int[n];
//		int numEdges = (int) ((Math.pow(n,2)) - n)/2;
//		Edge[] edges = new Edge[numEdges];
//		PriorityQueue q = new PriorityQueue(numEdges);
//		double time = 0;
//		double primSumL = 0;
//		
//		//create array of edges
//		/*int nextEdge = 0;
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < i; j++) {
//				edges[nextEdge] = new Edge(g[i][j], i, j);
//				nextEdge++;
//			}
//		}*/
//		
//		//select random vertice to start at
//		int randVertice = (int) (Math.random() * (n-1));
//		//tree = new DisjointSet(randVertice);
//		tree[randVertice]++; //mark this node as added to our MST
//		//initialise priority queue with edges connected to the randomly selected starting vertice
//		/*for(int i = 0; i < numEdges; i++) {
//			if((randVertice == edges[i].src) || (randVertice == edges[i].dest)) {
//				q.insert(edges[i]);
//			}
//		}*/
//		for(int j = 0; j < n; j++) {
//			if(tree[j] == 0) {
//				Edge connectedEdge = new Edge(g[randVertice][j], randVertice, j);
//				q.insert(connectedEdge);
//			}
//		}
//		
//		while(!q.isEmpty()) {
//			Edge e = q.extractMin();
//			int vertice = -1;
//			if((tree[e.src] == 1) && !(tree[e.dest] == 1)) {
//				vertice = e.dest;
//			}else if((tree[e.dest] == 1) && !(tree[e.src] == 1)){
//				vertice = e.src;
//			}
//			if(vertice != -1) {
//				primSumL += e.weight;
//				
//				if(n == 1000) {
//					System.out.println("Prim edge: "+e.weight);
//				}
//				
//				tree[vertice]++; //mark this vertice as added to our MST
//				/*for(int i = 0; i < numEdges; i++) {
//					if((vertice == edges[i].src) && !(tree.has(edges[i].dest))) {
//						q.insert(edges[i]);
//					}else if((vertice == edges[i].dest) && !(tree.has(edges[i].src))) {
//						q.insert(edges[i]);
//					}
//				}*/
//				for(int j = 0; j < n; j++) {
//					if(tree[j] == 0) {
//						Edge connectedEdge = new Edge(g[randVertice][j], randVertice, j);
//						q.insert(connectedEdge);
//					}
//				}
//			}
//		}
//		
//		//primAveL = primAveL / (n-1); //n-1 == number of edges in MST
//		/*System.out.print("Nodes in Prim MST: ");
//		for(int i = 0; i < n; i++) {
//			System.out.print(tree[i]+", ");
//		}*/
//		return primSumL;
//	}
	
	private static double prim(double[][] g) {
		int n = g[0].length;
		double[][] mst = new double[n][n];
		PriorityQueue q = new PriorityQueue(n);
		int[] tree = new int[n];
		double primSumL = 0;
		
		for(int i = 0; i < n; i++) {
			PQVert v = new PQVert(i);
			q.insert(v);
		}
		
		int randVertice = (int) (Math.random() * (n-1));
		System.out.println("Rand: "+randVertice);
		//q.queue[randVertice+1].key = 0;
		System.out.println("Vert: "+randVertice+", pos: "+q.get(randVertice)+" key: "+q.queue[q.get(randVertice)].key);
		q.decreaseKey(q.get(randVertice), 0);
		System.out.println("Vert: "+randVertice+", pos: "+q.get(randVertice)+" key: "+q.queue[q.get(randVertice)].key);
		while(!q.isEmpty()) {
			PQVert v = q.extractMin();
			tree[v.id] = 1;
			
			//add edge to mst
			primSumL += v.key;
			
			for(int j = 0; j < n; j++) {
				if((tree[j] != 1) && (g[v.id][j] != 0)) {
					int vert = q.get(j);
					if(g[v.id][j] < q.queue[vert].key) {
						q.queue[vert].pi = v;
						q.decreaseKey(vert, g[v.id][j]);
					}
				}
			}
		}
		return primSumL;
	}
	
	private static Edge[] quickSort(Edge[] edges, int low, int high) {
		int i = low;
		int j = high;
		double pivot = edges[i + ((j-i)/2)].weight;
		
		while(i <= j) {
			while(edges[i].weight < pivot) {
				i++;
			}
			while(edges[j].weight > pivot) {
				j--;
			}
			
			//swap values
			if(i <= j) {
				Edge temp = edges[i];
				edges[i] = edges[j];
				edges[j] = temp;
				
				i++;
				j--;
			}
		}
		
		if(low < j) edges = quickSort(edges, low, j);
		if(i < high) edges = quickSort(edges, i, high);
		
		return edges;
	}
	
	/* Directed Forest of DisjointSet Operations
	 * ----------------------------------------
	 -------------------------------------------*/
	public static void union(DisjointSet tree1, DisjointSet tree2) {
		if(tree1.rank > tree2.rank) {
			link(tree2, tree1);
		}else if(tree2.rank > tree1.rank) {
			link(tree1, tree2);
		}else {
			tree1.rank++;
			link(tree2, tree1);
		}
	}
	
	public static DisjointSet findSet(DisjointSet v) {
		while(v.parent != v) {
			v = v.parent;
		}
		return v;
	}
	
	public static void link(DisjointSet tree1, DisjointSet tree2) {
		tree1.parent = tree2;
	}
	/*--------------------------------------------*/
}

//Priority queue vertice
class PQVert {
	int id;
	PQVert pi;
	double key;
	
	public PQVert(int v) {
		id = v;
		pi = null;
		key = Double.MAX_VALUE;
	}
}

/*Priority Queue */
class PriorityQueue {
	PQVert[] queue;
	int size;
	
	//constructor
	public PriorityQueue(int n) {
		queue = new PQVert[n+1];
		size = 0;
	}
	
	//returns the index of a node with a given 'id' value
	public int get(int v) {
		for(int i = 1; i <= size; i++) {
			if(queue[i].id == v) {
				return i;
			}
		}
		System.out.println("Can't find "+v);
		return -1;
	}
	
	//inserts a single node into the queue
	public void insert(PQVert v) {
		queue[size + 1] = v;
		size++;
		fixQueueInsert();
	}
	
	//extracts node with the minimum key in the queue (in index 1)
	public PQVert extractMin() {
		System.out.println("Extracted "+queue[1].id);
		PQVert v = queue[1];
		queue[1] = queue[size];
		
		queue[size] = null;
		
		fixQueueExtract();
		//size--;
		System.out.print("Queue after fix: ");
		for(int i = 1; i < size; i++) {
			System.out.print(queue[i].id+", ");
		}
		System.out.println("");
		return v;
	}
	
	//decreases the key of a given vertice in the priority queue, and swaps nodes in the queue to maintain the fundamental properties of the min-heap priority queue
	public void decreaseKey(int v, double key) {
		queue[v].key = key;
		while((v > 1) && (queue[(int) (v / 2)].key > queue[v].key)) {
			System.out.println("Swapping "+(int) (v / 2)+"  with  "+v);
			PQVert temp = queue[(int) (v / 2)];
			queue[(int) (v / 2)] = queue[v];
			queue[v] = temp;
			v = (int) (v / 2);
		}
	}
	
	//returns true if the queue is empty
	public boolean isEmpty() {
		return (size == 0);
	}
	
	//adjusts queue to maintain properties after an insert operation
	private void fixQueueInsert() {
		int newVertIndex = size;
		PQVert v = queue[newVertIndex];
		
		if (newVertIndex == 1) return;
		
		int parent = (int) (newVertIndex / 2);
		while (queue[parent].key > v.key) {
			PQVert temp = queue[parent];
			queue[parent] = v;
			queue[newVertIndex] = temp;
			
			newVertIndex = parent;
			parent = (int) (newVertIndex / 2);
			if(newVertIndex == 1) return;
		}
	}
	
	//adjusts queue to maintain properties after an extract operation
	private void fixQueueExtract() {
		int parent = 1;
		int left = parent * 2;
		int right = (parent * 2) + 1;
		boolean notFixed = true;
		
		while(notFixed) {
			System.out.println("parent: "+parent+" left: "+left+" right: "+right);
			System.out.print("     parent: "+queue[parent].key);
			System.out.print("     left: "+queue[left].key);
			System.out.println("     right: "+queue[right].key);
			notFixed = false;
			if(queue[left] == null) {
				return; //the bottom of the priority queue has been reached;
			}else if((queue[right] == null) && (queue[left] != null)){
				if(queue[parent].key > queue[left].key) {
					PQVert temp = queue[parent];
					queue[parent] = queue[left];
					queue[left] = temp;
					parent = left;
					return; //the bottom of the priority queue has been reached
				}
			}else if(queue[left].key <= queue[right].key) {
				
				if(queue[parent].key > queue[left].key) {
					notFixed = true;
					PQVert temp = queue[parent];
					queue[parent] = queue[left];
					queue[left] = temp;
					parent = left;
				}
			}else {
				if(queue[parent].key > queue[right].key) {
					notFixed = true;
					PQVert temp = queue[parent];
					queue[parent] = queue[right];
					queue[right] = temp;
					parent = right;
				}
			}
			
			left = parent * 2;
			right = (parent * 2) + 1;
		}
	}
}

//Disjoint Set
class DisjointSet {
	DisjointSet parent;
	int rank;
	int value;
	
	public DisjointSet(int v) {
		parent = this;
		value = v;
		rank = 0;
	}
}

//Edge
class Edge {
	double weight;
	int src;
	int dest;

	//constructor
	public Edge(double w, int s, int d) {
		weight = w;
		src = s;
		dest = d;
	}
}