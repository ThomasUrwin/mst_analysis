package comp3600_ass_2;

import java.util.ArrayList;
import java.util.List;

public class mst {

	public static void main(String[] args) {
		System.out.println("Running test...");
		test();
	}
	
	private static void test() {
		Vertice v1 = new Vertice(1);
		Vertice v2 = new Vertice(2);
		System.out.println("V1 id: " + v1.id);
		System.out.println("V2 id: " + v2.id);
		
		Edge e = new Edge((float)0.5, v1, v2);
		
		System.out.println("E weight: " + e.weight + ", E source: " + e.src.id + ", E destination: " + e.dest.id);
	}
	
}

class Vertice {
	int id;
	
	//constructor
	public Vertice(int i) {
		id = i;
	}
}

class Edge {
	float weight;
	Vertice src;
	Vertice dest;
	
	//constructor
	public Edge(float w, Vertice s, Vertice d) {
		weight = w;
		src = s;
		dest = d;
	}
}

class Graph {
	List<Vertice> v;
	List<Edge> e;
	
	public Graph() {
		v = new ArrayList<Vertice>();
		e = new ArrayList<Edge>();
	}
	
	public void addV(Vertice vert) {
		v.add(vert);
	}
	
	public void removeV(int id) {
		
	}
	
	public void addE(Edge edge) {
		e.add(edge);
	}
	
	public void removeE(Edge edge) {
		
	}
}