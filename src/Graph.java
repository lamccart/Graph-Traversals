import sun.jvm.hotspot.oops.Array;

import java.util.*;

public class Graph {

    HashMap<Integer, Vertex> vertexMap;

    /**
     * Constructor for Graph
     */
    public Graph() {

        vertexMap = new HashMap<>();

    }

    /**
     * Adds a vertex to the graph. Throws IllegalArgumentException if given vertex
     * already exist in the graph.
     *
     * @param v vertex to be added to the graph
     * @throws IllegalArgumentException if two vertices with the same name are added.
     */
    public void addVertex(Vertex v) throws IllegalArgumentException {

        int hashKey = v.hashCode();

        if(vertexMap.containsKey(hashKey)){
            throw new IllegalArgumentException();
        }else{
            vertexMap.put(hashKey, v);
        }

    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {

        return vertexMap.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {

        int hashKey = name.hashCode();
        return vertexMap.get(hashKey);
    }

    /**
     * Adds a directed edge from vertex u to vertex v, Throws IllegalArgumentException if one of
     * the vertex does not exist
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight weight of the edge between vertex u and v
     * @throws IllegalArgumentException if one of the vertex does not exist
     */
    public void addEdge(String nameU, String nameV, Double weight) throws IllegalArgumentException {

        int hashKeyU = nameU.hashCode();
        int hashKeyV = nameV.hashCode();

        if(!vertexMap.containsKey(hashKeyU) || !vertexMap.containsKey(hashKeyV)){
            throw new IllegalArgumentException();
        }else{
            Vertex U = vertexMap.get(hashKeyU);
            Vertex V = vertexMap.get(hashKeyV);
            Edge newEdge = new Edge(U, V, weight);
            U.setEdge(newEdge);
        }

    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a directed
     * edge from u to v, then a directed edge from v to u
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight  weight of the edge between vertex u and v
     */
    public void addUndirectedEdge(String nameU, String nameV, double weight) {

        addEdge(nameU, nameV, weight);
        addEdge(nameV, nameU, weight);
    }

    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux (double) x coordinate of point u
     * @param uy (double) y coordinate of point u
     * @param vx (double) x coordinate of point v
     * @param vy (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {

        return Math.sqrt(Math.pow(vx - ux, 2) + Math.pow(vy - uy, 2));
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    public void computeAllEuclideanDistances() {

        //Loop through each vertex in map
        for(Vertex v : getVertices()){
            //Loop through each vertex edge has
            for(Edge e : v.getEdges()){
                //Compute euclidean distance between source and target
                double euclidDist = computeEuclideanDistance(e.getSource().getX(), e.getSource().getY(),
                                                            e.getTarget().getX(), e.getTarget().getY());
                //Set distance of current edge to euclidean distance
                e.setDistance(euclidDist);
            }
        }

    }

    /**
     * Helper method to reset all the vertices before doing graph traversal algorithms
     */
    private void resetAllVertices() {

        for(Vertex v : getVertices()){
            v.setVisited(false);
            v.setPredEdge(null);
            v.setDistance(Double.MAX_VALUE);
        }

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using DFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void DFS(String s, String t) {
        Stack<Vertex> fronteir = new Stack<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);

        fronteir.push(startV);
        startV.setVisited(true);
        while(!fronteir.empty()){
            Vertex currentV = fronteir.pop();
            for(Edge e : currentV.getEdges()){
                if(!e.getTarget().wasVisited()){
                    e.getTarget().setPredEdge(e);
                    fronteir.push(e.getTarget());
                }
            }
            currentV.setVisited(true);
            if(currentV == endV){
                break;
            }

        }

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using BFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void BFS(String s, String t) {
        LinkedList<Vertex> fronteir = new LinkedList<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);

        fronteir.push(startV);
        startV.setVisited(true);

        while(!fronteir.isEmpty()){
            Vertex currentV = fronteir.poll();
            if(currentV == endV){
                break;
            }
            for(Edge e : currentV.getEdges()){
                if(!e.getTarget().wasVisited()){
                    e.getTarget().setVisited(true);
                    e.getTarget().setPredEdge(e);
                    fronteir.add(e.getTarget());
                }
            }
        }

    }

    /**
     * Helper class for Dijkstra and A*, used in priority queue
     */
    private class CostVertex implements Comparable<CostVertex> {
        double cost;
        Vertex vertex;

        public CostVertex(double cost, Vertex vertex) {
            this.cost = cost;
            this.vertex = vertex;
        }

        public int compareTo(CostVertex o) {
            return Double.compare(cost, o.cost);
        }
    }

    /**
     * Find the shortest path from vertex with name s to vertex with name t, using Dijkstra
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void Dijkstra(String s, String t) {
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);
        startV.setDistance(0);
        startV.setPredEdge(null);
        queue.add(new CostVertex(startV.getDistance(), startV));

        while(!queue.isEmpty()){
            CostVertex currentV = queue.poll();

            if(currentV.vertex == endV){
                break;
            }

            for(Edge e : currentV.vertex.getEdges()){
                double newDistance = currentV.cost + e.getDistance();
                if(!e.getTarget().wasVisited() || newDistance < e.getTarget().getDistance()){
                    e.getTarget().setDistance(newDistance);
                    queue.add(new CostVertex(newDistance, e.getTarget()));
                    e.getTarget().setPredEdge(e);
                    e.getTarget().setVisited(true);
                }
            }
        }

    }

    /**
     * Helper method to calculate the h value in A*
     *
     * @param cur the current vertex being explored
     * @param goal the goal vertex to reach
     * @return the h value of cur and goal vertices
     */
    private double hValue(String cur, String goal) {

        Vertex curV = vertexMap.get(cur.hashCode());
        Vertex goalV = vertexMap.get(goal.hashCode());

        return computeEuclideanDistance(curV.getX(), curV.getY(), goalV.getX(), goalV.getY());
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using A*
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void AStar(String s, String t) {

        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);
        startV.setDistance(0);
        startV.setPredEdge(null);
        queue.add(new CostVertex(startV.getDistance(), startV));

        while(!queue.isEmpty()){
            CostVertex currentV = queue.poll();

            if(currentV.vertex == endV){
                break;
            }

            for(Edge e : currentV.vertex.getEdges()){
                double newDistance = currentV.vertex.getDistance() + e.getDistance();
                if(!e.getTarget().wasVisited() || newDistance < e.getTarget().getDistance()){
                    e.getTarget().setDistance(newDistance);
                    queue.add(new CostVertex(newDistance + hValue(t, e.getTarget().getName()), e.getTarget()));
                    e.getTarget().setPredEdge(e);
                    e.getTarget().setVisited(true);
                }
            }
        }

    }

    /**
     * Returns a list of edges for a path from city s to city t.
     *
     * @param s starting city name
     * @param t ending city name
     * @return list of edges from s to t
     */
    public List<Edge> getPath(String s, String t) {
        Vertex startV = vertexMap.get(s.hashCode());
        Vertex endV = vertexMap.get(t.hashCode());
        LinkedList<Edge> path = new LinkedList<>();
        if(endV.wasVisited()){
            Edge currentEdge = endV.getPredEdge();
            while(currentEdge != startV.getPredEdge()){
                path.add(currentEdge);
                currentEdge = currentEdge.getSource().getPredEdge();
            }
        }
        resetAllVertices();
        return path;
    }

}
