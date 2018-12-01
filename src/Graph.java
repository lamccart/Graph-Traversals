/**
 * Name: Liam McCarthy
 * PID: A14029718
 * Since: 11/25/2018
 */
import java.util.*;

public class Graph {

    HashMap<Integer, Vertex> vertexMap; //hashmap to hold all the graphs vertices

    /**
     * Constructor for Graph
     */
    public Graph() {
        //Initialize hashmap
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
        //Check if graph contains a vertex
        if(vertexMap.containsKey(hashKey)){
            throw new IllegalArgumentException();
        }else{
            //if not, insert it into the hashmap
            vertexMap.put(hashKey, v);
        }

    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {
        //Return the vertices of the hashmap
        return vertexMap.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {
        //return a vertex of a certain name
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
        //Check that both vertex names are in the hashmap
        if(!vertexMap.containsKey(hashKeyU) || !vertexMap.containsKey(hashKeyV)){
            throw new IllegalArgumentException();
        }else{
            //Get the vertices from hashmap
            Vertex U = vertexMap.get(hashKeyU);
            Vertex V = vertexMap.get(hashKeyV);
            //Make a new edge object
            Edge newEdge = new Edge(U, V, weight);
            //Set the directed edge on U
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
        //Set two directed edges, one from U to V and one from V to U
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
        //Compute the euclid distance between 2 points (x,y)
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
        //Loop through each vertex
        for(Vertex v : getVertices()){
            //Reset instance variables that are set in a graph transversal
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
        //Initialize fronteir and start and end vertices
        Stack<Vertex> fronteir = new Stack<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);

        //Push start vertex onto the stack and "visit" it
        fronteir.push(startV);
        startV.setVisited(true);
        //While there are vertices in the fronteir
        while(!fronteir.empty()){
            //Pop the top of stack
            Vertex currentV = fronteir.pop();
            //Loop through each edge of current vertex
            for(Edge e : currentV.getEdges()){
                //check if the vertex on other end of edge has been visited
                if(!e.getTarget().wasVisited()){
                    //if not visited, set the predecessor edge and push the target onto the stack
                    e.getTarget().setPredEdge(e);
                    fronteir.push(e.getTarget());
                }
            }
            //Set current vertex as visted
            currentV.setVisited(true);
            //See if current vertex is end vertex
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
        //Initialize fronteir and start and end vertex
        LinkedList<Vertex> fronteir = new LinkedList<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);

        //Push start vertex into queue and set as visited
        fronteir.push(startV);
        startV.setVisited(true);

        //Check if fronteir has vertices
        while(!fronteir.isEmpty()){
            //Take the first vertex in the queue
            Vertex currentV = fronteir.poll();
            //Check if it's the end vertex
            if(currentV == endV){
                break;
            }
            //Loop through edges connected to current vertex
            for(Edge e : currentV.getEdges()){
                //Check if edges target has been visited
                if(!e.getTarget().wasVisited()){
                    //If not then set it as visited
                    e.getTarget().setVisited(true);
                    //Set predecessor edge as current edge
                    e.getTarget().setPredEdge(e);
                    //Add target to the fronteir
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
        //Initialize priority queue and set start and end vertices
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);
        //Set start vertex cost to 0
        startV.setDistance(0);
        //Set it's predecessor edge to null
        startV.setPredEdge(null);
        //Add a new CostVertex to the priority queue for start vertex
        queue.add(new CostVertex(startV.getDistance(), startV));

        //Check if queue is empty
        while(!queue.isEmpty()){
            //Get top of heap
            CostVertex currentV = queue.poll();

            //Check it's not the end vertex
            if(currentV.vertex == endV){
                break;
            }
            //Loop through edges of current vertex
            for(Edge e : currentV.vertex.getEdges()){
                //Set the cost
                double newDistance = currentV.cost + e.getDistance();
                //Check if target was visited and if cost is less than target's cost
                if(!e.getTarget().wasVisited() || newDistance < e.getTarget().getDistance()){
                    //Set the targets cost the the new cost
                    e.getTarget().setDistance(newDistance);
                    //Add new cost vertex with target and new cost
                    queue.add(new CostVertex(newDistance, e.getTarget()));
                    //Set visited and precessor edge of target
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
        //Get current and goal vertices
        Vertex curV = vertexMap.get(cur.hashCode());
        Vertex goalV = vertexMap.get(goal.hashCode());
        //Compute their Euclidean distance
        return computeEuclideanDistance(curV.getX(), curV.getY(), goalV.getX(), goalV.getY());
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using A*
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void AStar(String s, String t) {
        //Initialize priority queue and set start and end vertices
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        Vertex startV = getVertex(s);
        Vertex endV = getVertex(t);
        //Set start vertex cost to 0
        startV.setDistance(0);
        //Set it's predecessor edge to null
        startV.setPredEdge(null);
        //Add a new CostVertex to the priority queue for start vertex
        queue.add(new CostVertex(startV.getDistance(), startV));

        //Check if queue is empty
        while(!queue.isEmpty()){
            //Get top of heap
            CostVertex currentV = queue.poll();

            //Check it's not the end vertex
            if(currentV.vertex == endV){
                break;
            }

            //Loop through edges of current vertex
            for(Edge e : currentV.vertex.getEdges()){
                //Set the cost
                double newDistance = currentV.vertex.getDistance() + e.getDistance();
                //Check if target was visited and if cost is less than target's cost
                if(!e.getTarget().wasVisited() || newDistance < e.getTarget().getDistance()){
                    //Set the targets cost the the new cost
                    e.getTarget().setDistance(newDistance);
                    //Add new cost vertex with target and new cost + the hValue
                    queue.add(new CostVertex(newDistance + hValue(t, e.getTarget().getName()), e.getTarget()));
                    //Set visited and predecessor edge of target
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
        //Get the start and end vertices
        Vertex startV = vertexMap.get(s.hashCode());
        Vertex endV = vertexMap.get(t.hashCode());
        //Initialize path
        LinkedList<Edge> path = new LinkedList<>();
        //Check if the end vertex was visited
        if(endV.wasVisited()){
            //Get the predecessor edge to the end vertex
            Edge currentEdge = endV.getPredEdge();
            //Check if we have not reached the start
            while(currentEdge != startV.getPredEdge()){
                //add the current edge to path
                path.add(currentEdge);
                //set current edge to the predecessor edge of the source
                currentEdge = currentEdge.getSource().getPredEdge();
            }
        }
        //Reset all instance variables and return path
        resetAllVertices();
        return path;
    }

}
