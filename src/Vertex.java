/**
 * Name: Liam McCarthy
 * PID: A14029718
 * Since: 12/1/2018
 */
import java.util.LinkedList;

public class Vertex {

    private final String name; // the name of this vertex
    private final int x; // the x coordinates of this vertex on map
    private final int y; // the y coordinates of this vertex on map
    private boolean visited; // checks if vertex was visited
    private Edge predEdge; // sets the edge visted before this vertex
    private double distance; // sets the cost of the vertex
    private LinkedList<Edge> edges; //edges the vertex is connected to

    public Vertex(String name, int x, int y) {
        //Initialize instance variables
        this.name = name;
        this.x = x;
        this.y = y;
        edges = new LinkedList<>();
        visited = false;
        predEdge = null;
        distance = Double.MAX_VALUE;
    }

    public String getName() {
        return name;
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public LinkedList<Edge> getEdges(){
        return edges;
    }

    public void setEdge(Edge e){
        edges.addLast(e);
    }

    public boolean wasVisited(){
        return visited;
    }

    public void setVisited(boolean b){
        visited = b;
    }

    public double getDistance(){
        return distance;
    }

    public void setDistance(double d){
        distance = d;
    }

    public Edge getPredEdge(){
        return predEdge;
    }

    public void setPredEdge(Edge e){
        predEdge = e;
    }

    @Override
    public int hashCode() {
        // we assume that each vertex has a unique name
        return name.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null) {
            return false;
        }
        if (!(o instanceof Vertex)) {
            return false;
        }
        Vertex oVertex = (Vertex) o;

        return name.equals(oVertex.name) && x == oVertex.x && y == oVertex.y;
    }

    public String toString() {
        return name + " (" + x + ", " + y + ")";
    }

}