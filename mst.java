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
		double krusAve = 0;
		double krusL = 0;
		double primAve = 0;
		double primL= 0;
		double[][][] graphs = new double[20][n][n];
		
		//generate all 20 graphs, stored into an array
		for (int i = 0; i < 20; i++) {
			graphs[i] = initGraph(n);
		}
		System.out.println("Starting tests for n = "+n);
		for (int i = 0; i < 20; i++) {
			double[] krusVals;
			double[] primVals;
			long startTime = 0;
			long endTime = 0;
			
			//run Kruskal's and Prim's algorithms on the graph, returning average weight L(n) and time to complete
			krusVals = kruskal(graphs[i]);
			startTime = System.currentTimeMillis();
			primVals = prim(graphs[i]);
			endTime = System.currentTimeMillis();
			
			krusAve += krusVals[0];
			krusL += krusVals[1];
			primAve += (endTime - startTime);
			primL += primVals[1];
			
			System.out.println("Iteration done... "+n+": "+i);
		}
		
		//Calculate averages
		krusAve = krusAve / 20;
		krusL = krusL / 20;
		primAve = primAve / 20;
		primL = primL / 20;
		
		System.out.println("n = " + String.valueOf(n) + " --> KrusAve: " + krusAve + ", KrusL: " + krusL + ", PrimAve: " + primAve + ", PrimL: " + primL);
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
	private static double[] kruskal(double[][] g) {
		int n = g[0].length;
		double[][] mst = new double[n][n];
		List<DisjointSet> forest = new ArrayList<DisjointSet>();
		int numEdges = (int) ((Math.pow(n,2)) - n)/2;
		Edge[] edges = new Edge[numEdges];
		double time = 0;
		double krusL = 0;
		
		long startTime = System.currentTimeMillis();
		
		//create DisjointSet for each vertex
		for (int i = 0; i < n; i++) {
			DisjointSet t = new DisjointSet(i);
			forest.add(t);
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
			
			for(DisjointSet t : forest) {
				if(t.has(edges[i].src)) {
					treeSrc = t;
				}
				if(t.has(edges[i].dest)) {
					treeDest = t;
				}
				
				if((treeSrc != null) && (treeDest != null)) {
					if(treeSrc.getRep() != treeDest.getRep()) {
						//add edge to MST
						mst[edges[i].src][edges[i].dest] = edges[i].weight;
						mst[edges[i].dest][edges[i].src] = edges[i].weight;
						
						krusL += edges[i].weight;
						
						//union sets
						treeSrc = unionDisjSet(treeSrc, treeDest);
						forest.remove(treeDest);
						break;
					}else break;
				}
			}
		}
		long endTime = System.currentTimeMillis();

		time = endTime - startTime;
		krusL = krusL / (n-1); //n-1 == number of edges in MST
		
		double[] returnVals = new double[2];
		returnVals[0] = time;
		returnVals[1] = krusL;
		return returnVals;
	}
	
	
	//PRIM'S ALGORITHM
	private static double[] prim(double[][] g) {
		int n = g[0].length;
		double[][] mst = new double[n][n];
		DisjointSet tree = null;
		int numEdges = (int) ((Math.pow(n,2)) - n)/2;
		Edge[] edges = new Edge[numEdges];
		PriorityQueue q = new PriorityQueue(numEdges);
		double time = 0;
		double primL = 0;
		
		//create array of edges
		int nextEdge = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				edges[nextEdge] = new Edge(g[i][j], i, j);
				nextEdge++;
			}
		}
		
		//select random vertice to start at
		int randVertice = (int) (Math.random() * (n-1));
		tree = new DisjointSet(randVertice);
		//initialise priority queue with edges connected to the randomly selected starting vertice
		for(int i = 0; i < numEdges; i++) {
			if((randVertice == edges[i].src) || (randVertice == edges[i].dest)) {
				q.enqueue(edges[i]);
			}
		}
		
		while(!q.isEmpty()) {
			Edge e = q.dequeue();
			int vertice = -1;
			if((tree.has(e.src)) && !(tree.has(e.dest))) {
				vertice = e.dest;
			}else if((tree.has(e.dest)) && !(tree.has(e.src))){
				vertice = e.src;
			}
			if(vertice != -1) {
				primL += e.weight;
				tree.include(vertice);
				for(int i = 0; i < numEdges; i++) {
					if((vertice == edges[i].src) && !(tree.has(edges[i].dest))) {
						q.enqueue(edges[i]);
					}else if((vertice == edges[i].dest) && !(tree.has(edges[i].src))) {
						q.enqueue(edges[i]);
					}
				}
			}
		}
		
		primL = primL / (n-1); //n-1 == number of edges in MST
		
		double[] returnVals = new double[2];
		returnVals[0] = time;
		returnVals[1] = primL;
		return returnVals;
	}
	
	private static DisjointSet unionDisjSet(DisjointSet t1, DisjointSet t2) {
		DisjointSet newSet;
		DisjointSet childSet;
		if (t1.size >= t2.size) {
			newSet = t1;
			childSet = t2;
		}else {
			newSet = t2;
			childSet = t1;
		}	
		newSet.include(childSet);
		
		return newSet;
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
}

class PriorityQueue {
	
}
/*class PriorityQueue {
	Edge[] queue;
	int front;
	int rear;
	
	public PriorityQueue(int n) {
		queue = new Edge[n];
		front = 0;
		rear = 0;
	}
	
	private void reorder() {
		Edge e = queue[rear];
		for(int i = front; i < rear; i++) {
			if(queue[i].weight > e.weight) {
				if((front > 0) && (i == front)) {
					queue[front-1] = e;
					queue[rear] = null;
					front--;
					rear--;
				}else {
					for(int j = rear-1; j >= i; j--) {
						if(j != -1) {
							queue[j+1] = queue[j];
						}else {
							break;//there is only 1 edge in the queue, no reordering is needed
						}
						
					}
					queue[i] = e;
				}
				break;
			}
		}
	}
	
	public void enqueue(Edge e) {
		if(queue[rear] != null) {
			queue[rear+1] = e;
		}else {
			queue[rear] = e;
		}
		if(queue[rear+1] != null) rear++;
		reorder();
	}
	
	public Edge dequeue() {
		Edge e = queue[front];
		queue[front] = null;
		if(front < rear) front++;
		return e;
	}
	
	public boolean isEmpty() {
		if((front == rear) && (queue[front] == null)) {
			return true;
		}else {
			return false;
		}
	}
}*/

class DisjointSet {
	int rep;
	Set<Integer> vertices;
	int size;
	
	//constructor
	public DisjointSet(int v) {
		vertices = new HashSet<Integer>();
		rep = v;
		vertices.add(v);
		size = 1;
	}
	
	//checks if a given vertice is in the set
	public boolean has(int v) {
		return vertices.contains(v);
	}
	
	//adds a new vertices to the set
	public void include(int v) {
		vertices.add(v);
		size += 1;
	}
	
	public void include(DisjointSet t) {
		for (int v : t.vertices) {
			this.include(v);
		}
	}
	
	//returns the representative of the set
	public int getRep() {
		return rep;
	}
}

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


//public class mst {
//
//	public static void main(String[] args) {
//		System.out.println("Running test...");
//		test();
//	}
//	
//	private static void test() {
//		Vertice v1 = new Vertice(1);
//		Vertice v2 = new Vertice(2);
//		System.out.println("V1 id: " + v1.id);
//		System.out.println("V2 id: " + v2.id);
//		
//		Edge e = new Edge((float)0.5, v1, v2);
//		
//		System.out.println("E weight: " + e.weight + ", E source: " + e.src.id + ", E destination: " + e.dest.id);
//	}
//	
//}
//
//class Vertice {
//	int id;
//	
//	//constructor
//	public Vertice(int i) {
//		id = i;
//	}
//}
//
//class Edge {
//	float weight;
//	Vertice src;
//	Vertice dest;
//	
//	//constructor
//	public Edge(float w, Vertice s, Vertice d) {
//		weight = w;
//		src = s;
//		dest = d;
//	}
//}
//
//class Graph {
//	List<Vertice> v;
//	List<Edge> e;
//	
//	public Graph() {
//		v = new ArrayList<Vertice>();
//		e = new ArrayList<Edge>();
//	}
//	
//	public void addV(Vertice vert) {
//		v.add(vert);
//	}
//	
//	public void removeV(int id) {
//		
//	}
//	
//	public void addE(Edge edge) {
//		e.add(edge);
//	}
//	
//	public void removeE(Edge edge) {
//		
//	}
//}