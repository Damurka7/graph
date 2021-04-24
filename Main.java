import java.util.*;
import java.util.stream.Collectors;

import static java.lang.Math.min;

public class Main {
    public static void main(String[] args) {
        //inputFirstSecondTask();
        inputThirdTask();
    }

    /**
     * input for the first two tasks
     * scanning from input while there is something to insert
     */
    public static void inputFirstSecondTask() {
        AdjacencyMatrixGraph graph = new AdjacencyMatrixGraph();
        int n;
        Scanner sc = new Scanner(System.in);
//        n = sc.nextInt();
        while (sc.hasNext()) {
            String s = sc.nextLine();
            String new_str[] = s.split(" ");
            String sub = new_str[0];
            switch (sub) {
                case "ADD_VERTEX":
                    graph.addVertex(new_str[1]);
                    //      graph.printAdjacencyMatrixGraph();
                    break;
                case "REMOVE_VERTEX":
                    graph.removeVertex(graph.findVertex(new_str[1]));
                    //       graph.printAdjacencyMatrixGraph();
                    break;
                case "ADD_EDGE":
                    //AdjacencyMatrixGraph.Edge e = new AdjacencyMatrixGraph.Edge(graph.findVertex(new_str[1]), graph.findVertex(new_str[2]), Integer.getInteger(new_str[3]));
                    graph.addEdge(graph.findVertex(new_str[1]), graph.findVertex(new_str[2]), Integer.parseInt(new_str[3]));
                    //  graph.printAdjacencyMatrixGraph();
                    break;
                case "REMOVE_EDGE":
                    AdjacencyMatrixGraph.Edge edge = new AdjacencyMatrixGraph.Edge(graph.findVertex(new_str[1]), graph.findVertex(new_str[2]));
                    graph.removeEdge(edge);
                    //    graph.printAdjacencyMatrixGraph();
                    break;
                case "HAS_EDGE":
                    // System.out.println(graph.hasEdge(graph.findVertex(new_str[1]), graph.findVertex(new_str[2])));
                    if (graph.hasEdge(graph.findVertex(new_str[1]), graph.findVertex(new_str[2]))) {
                        System.out.println("TRUE");
                    } else {
                        System.out.println("FALSE");
                    }
                    // graph.printAdjacencyMatrixGraph();
                    break;
                case "TRANSPOSE":
                    graph.trnspose();
                    // graph.printAdjacencyMatrixGraph();

                    break;
                case "IS_ACYCLIC":
//                    if (graph.isAcyclic()) {
//                        System.out.println("ACYCLIC");
//                    } else {
//                        List<AdjacencyMatrixGraph.Vertex> list = graph.findCycle();
//                        System.out.print(graph.getSum() + " ");
//                        for (int i = 0; i < list.size()-1; i++) {
//                            AdjacencyMatrixGraph.Vertex v = list.get(i);
//                            System.out.print(v.getName() + " ");
//                        }
//                        AdjacencyMatrixGraph.Vertex v = list.get(list.size()-1);
//                        System.out.println(v.getName() + " ");
//                    }

                    if (!graph.isAcyclic()) {
                        List<AdjacencyMatrixGraph.Vertex<String>> cycle = graph.findCycle();
                        long cycleWeight = 0;
                        for (int i = 0; i < cycle.size() - 1; i++) {
                            cycleWeight += (Integer)graph.findEdge(cycle.get(i), cycle.get(i + 1)).getWeight();
                        }
                        cycleWeight += (Integer)graph.findEdge(cycle.get(cycle.size() - 1), cycle.get(0)).getWeight();
                        System.out.print(cycleWeight + " ");
                        for (AdjacencyMatrixGraph.Vertex v: cycle){
                            System.out.print(v.getName() + " ");
                        }
                        System.out.println();
                    } else {
                        System.out.println("ACYCLIC");
                    }

                    //graph.printAdjacencyMatrixGraph();

                    break;
//                case "D":
//                    graph.Dijkstra(graph.findVertex(new_str[1]), graph.findVertex(new_str[2]));
//                   // graph.printAdjacencyMatrixGraph();
//                    //graph.pr();
//                    break;

            }
        }
    }

    /**
     * input for the third task
     */
    public static void inputThirdTask() {
        int n, m, w;
        AdjacencyMatrixGraph<Integer, Integer> graph = new AdjacencyMatrixGraph<>();
        Scanner sc = new Scanner(System.in);
        n = sc.nextInt();
        m = sc.nextInt();

        for (int i = 0; i < n; i++) {
            graph.addVertex(i + 1);
        }

        for (int i = 0; i < m; i++) {
            int begV = sc.nextInt();
            int destV = sc.nextInt();
            int weight = sc.nextInt();
            int band = sc.nextInt();
            graph.addEdge(graph.findVertex(begV), graph.findVertex(destV), weight, band);

        }
        n = sc.nextInt();
        m = sc.nextInt();
        w = sc.nextInt();

        graph.Dijkstra(graph.findVertex(n), graph.findVertex(m), w);

        // graph.printAdjacencyMatrixGraph();
    }
}

interface GraphADT<V extends Comparable<V>, E extends Comparable<E>> {
    AdjacencyMatrixGraph.Vertex<V> addVertex(V value);

    void removeVertex(AdjacencyMatrixGraph.Vertex<V> vertex);

    AdjacencyMatrixGraph.Edge<V, E> addEdge(AdjacencyMatrixGraph.Vertex<V> from, AdjacencyMatrixGraph.Vertex<V> to, E weight);

    void removeEdge(AdjacencyMatrixGraph.Edge<V, E> e);

    ArrayList<AdjacencyMatrixGraph.Edge<V, E>> edgesFrom(AdjacencyMatrixGraph.Vertex<V> v);

    ArrayList<AdjacencyMatrixGraph.Edge<V, E>> edgesTo(AdjacencyMatrixGraph.Vertex<V> v);

    AdjacencyMatrixGraph.Vertex<V> findVertex(V value);

    AdjacencyMatrixGraph.Edge<V, E> findEdge(AdjacencyMatrixGraph.Vertex<V> from, AdjacencyMatrixGraph.Vertex<V> to);

    boolean hasEdge(AdjacencyMatrixGraph.Vertex<V> u, AdjacencyMatrixGraph.Vertex<V> v);
}

class AdjacencyMatrixGraph<V extends Comparable<V>, E extends Comparable<E>> implements GraphADT<V, E> {
    private static final int INFINITY = 1000000000; //infinity - to find the smallest number (used in Dijkstra)

    //adjacency matrix which stores edge at each its cell
    public ArrayList<ArrayList<Edge<V, E>>> adjMatrix;

    //Each vertex has its own integer - position in adjacency matrix
    public HashMap<Vertex<V>, Integer> index;//to keep time complexity when finding edge
    public ArrayList<Vertex<V>> vertexList;

    public static class Vertex<T> implements Comparable<Vertex<T>> {
        private T name;
        private int distance; //used in Dijkstra to save distance from source to a current vertex
        public Vertex<T> parent; //used in Dijkstra to track the path from source to end vertex

        /**
         * constructor which creates a new vertex with name
         * @param name
         */
        public Vertex(T name) {
            this.name = name;
        }

        /**
         * @return a name of the vertex
         */
        public T getName() {
            return name;
        }

        /**
         * sets a distance
         * @param distance - how far this vertex from the source
         */
        public void setDistance(int distance) {
            this.distance = distance;
        }

        /**
         * @return a distance from source vertex
         */
        public int getDistance() {
            return distance;
        }

        /**
         * used in Priority Queue
         * each vertex is a pair <distance, Vertex>
         * @param o - vertex
         * @return - positive number if distance if this.vertex is greater tha a given, zero if they are equal, negative number otherwise
         */
        @Override
        public int compareTo(Vertex<T> o) {
            return this.distance - o.distance;
        }

        /**
         * compares two vertices
         * @param o - vertex
         * @return true if they are equal, false otherwise
         */
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Vertex<?> vertex = (Vertex<?>) o;

            return Objects.equals(name, vertex.name);
        }

        /**
         * @return a hashCode of given vertex by its name
         */
        @Override
        public int hashCode() {
            return this.name.hashCode();
        }
    }

    public static class Edge<V, E> {
        private E bandwidth;
        private E weight;
        private Vertex<V> destVertex;
        private Vertex<V> beginVertex;
        private boolean isNull;

        /**
         * constructor of an edge
         * @param begin - begin vertex
         * @param dest - destination vertex
         * @param w - weight of an edge
         * @param b - bandwidth of an edge
         */
        public Edge(Vertex<V> begin, Vertex<V> dest, E w, E b) {
            this.bandwidth = b;
            this.destVertex = dest;
            this.beginVertex = begin;
            this.weight = w;
        }

        /**
         * constructor of an edge
         * @param begin - begin vertex
         * @param dest - destination vertex
         * @param w - weight of an edge
         */
        public Edge(Vertex<V> begin, Vertex<V> dest, E w) {
            this.destVertex = dest;
            this.beginVertex = begin;
            this.weight = w;
        }

        /**
         * constructor of an edge
         * @param begin - begin vertex
         * @param dest - destination vertex
         */
        public Edge(Vertex<V> begin, Vertex<V> dest) {
            this.destVertex = dest;
            this.beginVertex = begin;
            this.isNull = false;
        }

        /**
         * constructor of an edge which creates a "null edge"
         */
        public Edge() {
            this.isNull = true;
        }

        /**
         * makes edge a "null edge"
         */
        public void remove() {
            this.isNull = true;
            this.destVertex = null;
            this.beginVertex = null;
            this.weight = null;
            this.bandwidth = null;

        }

        /**
         * @return weight of an edge
         */
        public E getWeight() {
            return weight;
        }

        /**
         * @return bandwidth of an edge
         */
        public E getBandwidth() {
            return bandwidth;
        }

        /**
         * @return destination vertex
         */
        public Vertex<V> getDestVertex() {
            return destVertex;
        }

        /**
         * @return a begin vertex
         */
        public Vertex<V> getBeginVertex() {
            return beginVertex;
        }
    }

    /**
     * constructor for graph
     */
    public AdjacencyMatrixGraph() {
        adjMatrix = new ArrayList<ArrayList<Edge<V, E>>>();
        index = new HashMap<Vertex<V>, Integer>();
        vertexList = new ArrayList<Vertex<V>>();
    }

    /**
     * creates an edge between given vertices and adds it to adjacency matrix
     * @param from - begin vertex
     * @param to - destination vertex
     * @param weight - weight of an edge
     * @param band - bandwidth
     * @return created edge
     */
    public Edge<V, E> addEdge(Vertex<V> from, Vertex<V> to, E weight, E band) {
        Edge<V, E> edge = new Edge<V, E>(from, to, weight, band);
        adjMatrix.get(index.get(from)).set(index.get(to), edge);
        return edge;
    }

    /**
     * adds new column and row to adjacency matrix
     * @param value - name of vertex
     * @return created vertex
     */
    @Override
    public Vertex<V> addVertex(V value) {
        if (!vertexList.contains(findVertex(value))) {
            Vertex<V> v = new Vertex<V>(value);
            this.index.put(v, this.vertexList.size());
            this.vertexList.add(v);//to maintain mapping of vertices and indexes
            ArrayList<Edge<V, E>> newRow = new ArrayList<>();

            Edge<V, E> edge = new Edge<>();
            for (int i = 0; i < adjMatrix.size(); i++) {
                newRow.add(edge);
            }
            adjMatrix.add(newRow);

            for (ArrayList<Edge<V, E>> matrix : adjMatrix) {
                Edge<V, E> e = new Edge<>();
                matrix.add(e);
            }
            return v;
        } else {
            return findVertex(value);
        }

    }

    /**
     * finds index of a given vertex
     * removes row and column with found index
     * removes this vertex from index and vertexList
     * changes the indices of vertices which are moved to the left in vertexList
     * @param vertex - given vertex
     */
    @Override
    public void removeVertex(Vertex<V> vertex) {
        int vertInd = index.get(vertex);
        for (ArrayList<Edge<V, E>> row : adjMatrix) {
            row.remove(vertInd);
        }
        adjMatrix.remove(vertInd);
        this.index.remove(vertex);
        for (int i = vertInd + 1; i < vertexList.size(); i++) {
            index.put(vertexList.get(i), i - 1);
            vertexList.set(i - 1, vertexList.get(i));
        }
        this.vertexList.remove(vertexList.size() - 1);
    }

    /**
     * finds the position of given vertices(constant time complexity) and adds edge at this cell
     * @param from - begin vertex
     * @param to - destination vertex
     * @param weight - weight of a given edge
     * @return created edge
     */
    @Override
    public Edge<V, E> addEdge(Vertex<V> from, Vertex<V> to, E weight) {
        Edge<V, E> edge = new Edge<V, E>(from, to, weight);
        adjMatrix.get(index.get(from)).set(index.get(to), edge);
        return edge;
    }

    /**
     * finds given edge (constant time complexity) and deletes it
     * @param e - given edge
     */
    @Override
    public void removeEdge(Edge<V, E> e) {
        adjMatrix.get(index.get(e.beginVertex)).get(index.get(e.destVertex)).remove();
    }

    /**
     * adds to the list row which corresponds to a given vertex
     * @param v - given vertex
     * @return a list of edges which are outgoes from a given vertex
     */
    @Override
    public ArrayList<Edge<V, E>> edgesFrom(Vertex<V> v) {
        ArrayList<Edge<V, E>> list = new ArrayList<>();
        ArrayList<Edge<V, E>> row = adjMatrix.get(index.get(v));
        for (Edge<V, E> edge : row) {
            if (!edge.isNull) {
                list.add(edge);
            }
        }
        return list;
    }

    /**
     * adds to the list column which corresponds to a given vertex
     * @param v - given vertex
     * @return a list of edges which are incomes to a given vertex
     */
    @Override
    public ArrayList<AdjacencyMatrixGraph.Edge<V, E>> edgesTo(Vertex<V> v) {
        ArrayList<Edge<V, E>> list = new ArrayList<>();
        for (ArrayList<Edge<V, E>> row : adjMatrix) {
            if (!row.get(index.get(v)).isNull) {
                list.add(row.get(index.get(v)));
            }
        }
        return list;
    }

    /**
     * returns a vertex with a given name
     * if there is no such vertex returns null
     * @param value - name of vertex
     * @return - vertex or null
     */
    @Override
    public Vertex<V> findVertex(V value) {
        for (Vertex<V> v : vertexList) {
            if (v.name.equals(value)) {
                return v;
            }
        }
        return null;
    }

    /**
     * finds edge in adjacency matrix
     * @param from - begin vertex
     * @param to - destination vertex
     * @return - found edge
     */
    @Override
    public Edge<V, E> findEdge(Vertex<V> from, Vertex<V> to) {
        return adjMatrix.get(index.get(from)).get(index.get(to));
    }

    /**
     * determine whether graph has edge between given vertices
     *
     * @param u given vertex
     * @param v given vertex
     * @return true if grapgh has edge, false - otherwise
     */
    @Override
    public boolean hasEdge(Vertex<V> u, Vertex<V> v) {
        if (vertexList.contains(u) && vertexList.contains(v)) {
            return !adjMatrix.get(index.get(u)).get(index.get(v)).isNull;
        } else {
            return false;
        }
    }

    /**
     * finds vertices which are adjacent to a given
     * stores them in list
     * @param vertex given vertex
     * @return list of adjacent vertices
     */
    public ArrayList<Vertex<V>> adjVertices(Vertex<V> vertex) {
        ArrayList<Edge<V, E>> adjEdges = edgesFrom(vertex);
        ArrayList<Vertex<V>> list = new ArrayList<>();
        for (Edge<V, E> edge : adjEdges)
            list.add(edge.destVertex);
        return list;
    }

    /**
     * prints the adjacency matrix
     */
    public void printAdjacencyMatrixGraph() {
        for (int i = 0; i < adjMatrix.size(); i++) {
            System.out.print("vertex name: " + vertexList.get(i).getName() + ": ");
            for (Edge<V, E> e : adjMatrix.get(i)) {
                if (!e.isNull) {
                    System.out.print(" destVertex_" + e.getDestVertex().getName() + " weight: " + e.getWeight() + " | ");
                } else {
                    System.out.print("|           null         |");
                }
            }
            System.out.print("\n");
        }
    }

    /**
     * is not used by user
     * method created to track correctness of Dijkstra algorithm
     */
    public void pr() {
        for (int i = 0; i < adjMatrix.size(); i++) {
            System.out.print("vertex name: " + vertexList.get(i).getName() + "; distance from A: " + vertexList.get(i).getDistance() + " |");
        }
    }

    /**
     * transposing adjacent matrix and reversing destination and beginning vertex for each edge
     */
    void trnspose() {
        ArrayList<ArrayList<Edge<V, E>>> nList = new ArrayList<>();
        for (int i = 0; i < adjMatrix.size(); i++) {
            ArrayList<Edge<V, E>> row = new ArrayList<>();
            for (int j = 0; j < adjMatrix.size(); j++) {
                Edge<V, E> edge = adjMatrix.get(j).get(i);
                if (!edge.isNull) {
                    Vertex<V> temp = edge.destVertex;
                    edge.destVertex = edge.beginVertex;
                    edge.beginVertex = temp;
                }
                row.add(edge);
            }
            nList.add(row);
        }
        this.adjMatrix = nList;
    }

    /**
     * checks whether the graph contains a cycle or not
     * @return false if graph has a cycle, true otherwise
     */
    public boolean isAcyclic() {
        return findCycle() == null;
    }

    /**
     * finds cycles with DFS
     * DFS is based on recursion (because we already did the stack version on labs)
     * @return list of cycle if there exists, otherwise empty list
     */
    public List<Vertex<V>> findCycle() {
        int[] p = new int[vertexList.size()];
        int[] colors = new int[vertexList.size()];
        List<Integer> list = null;
        for (Vertex<V> vertex : vertexList) {
            if (colors[index.get(vertex)] == 0) {
                list = new ArrayList<>();
                dfs(list, index.get(vertex), -1, p, colors);
                if (list.isEmpty()) {
                    list = null;
                } else {
                    break;
                }
            }
        }
        if (list == null) {
            return null;
        }
        ArrayList<Vertex<V>> cycle = new ArrayList<>();
        for (int index : list)
            cycle.add(vertexList.get(index));
        return cycle;
    }

    /**
     * DFS is based on recursion (because we already did the stack version on labs)
     * @param cycleList - list which contains a cycle
     * @param curVertex - vertex for which we are finding the cycle
     * @param prevVertex - previous vertex
     * @param arr - array which has previous vertex at currentVertex position
     * @param colors - array of colors for vertices
     */
    public void dfs(List<Integer> cycleList, int curVertex, int prevVertex, int[] arr, int[] colors) {
        colors[curVertex] = 1; // vertex is visited and is processing now
        arr[curVertex] = prevVertex;
        //ArrayList<Edge<V,E>> list = this.edgesFrom(vertexList.get(curVertex));
        for (Edge<V, E> edgeTo : this.edgesFrom(this.vertexList.get(curVertex))) {
            int to = this.index.get(edgeTo.getDestVertex());
            //int to = Integer.parseInt(String.valueOf(edge.destVertex));
            if (colors[to] == 0) {
                dfs(cycleList, to, curVertex, arr, colors);
                if (!cycleList.isEmpty()) {
                    break;
                }
            } else if (colors[to] == 1) {
                int current = curVertex;
                while (current != to) {
                    cycleList.add(current);
                    current = arr[current];
                }
                cycleList.add(current);
                List<Integer> newCycleList  = new ArrayList<>();
                //reversing
                for (int i = 0; i < cycleList.size(); i++) {
                    newCycleList.add(cycleList.get(cycleList.size()-i-1));
                }
                for (int i = 0; i < newCycleList.size(); i++) {
                    cycleList.set(i, newCycleList.get(i));
                }
                break;
            }
        }
        colors[curVertex] = 2;
    }

    /**
     * calculates the weight of the cycle
     * @return weight of the cycle
     */
    public int getSum() {
        List<Vertex<V>> list = findCycle();
        int sum = 0;

        if (!list.isEmpty()) {
            for (int i = 0; i < list.size() - 1; i++) {
                sum += (Integer) adjMatrix.get(index.get(list.get(i))).get(index.get(list.get(i + 1))).getWeight();
            }
            if (adjMatrix.get(index.get(list.get(0))).get(index.get(list.get(list.size() - 1))).isNull) {
                sum += (Integer) adjMatrix.get(index.get(list.get(0))).get(index.get(list.get(list.size() - 1))).getWeight();
            } else {
                sum += (Integer) adjMatrix.get(index.get(list.get(list.size() - 1))).get(index.get(list.get(0))).getWeight();
            }
        }

        return sum;
    }

    /**
     * extended Dijkstra algorithm which calculates shortest paths to each vertex from the source
     * for each vertex(each vertex has a list) method adds other vertices in its list which are on a path from source
     * @param source - starting vertex
     * @param end    - end vertex
     */
    public void Dijkstra(Vertex<V> source, Vertex<V> end, int w) {
        for (ArrayList<Edge<V, E>> row : adjMatrix) {
            for (Edge<V, E> edge : row) {
                if (!edge.isNull) {
                    if ((Integer) edge.bandwidth < w) {
                        edge.remove();
                    }
                }
            }
        }

        //printAdjacencyMatrixGraph();

        for (Vertex<V> vertex : vertexList) {
            if (vertex.equals(source)) {
                vertex.setDistance(0);
            } else {
                vertex.setDistance(INFINITY);
            }
        }

        PriorityQueue<Vertex<V>> queue = new PriorityQueue<>();
        queue.offer(source);
        while (!queue.isEmpty()) {
            Vertex<V> u = queue.remove();
            for (Edge<V, E> edge : edgesFrom(u)) {
                if (!edge.isNull) {
                    if (u.getDistance() + (Integer) edge.weight < edge.destVertex.distance) {
                        edge.destVertex.setDistance(u.distance + (Integer) edge.weight);
                        edge.destVertex.parent = u;
                        queue.offer(edge.destVertex);
                    }
                }
            }
        }

        if (end.distance != INFINITY) {
            ArrayList<Vertex<V>> reversedPath = new ArrayList<>();
            Vertex<V> curr = end;
            while (curr != null) {
                reversedPath.add(curr);
                curr = curr.parent;
            }
            ArrayList<Vertex<V>> path = new ArrayList<>();
            for (int i = reversedPath.size() - 1; i >= 0; i--) {
                path.add(reversedPath.get(i));
            }
            long length = 0;
            long bandwidth = INFINITY;
            for (int i = 0; i < path.size() - 1; i++) {
                Edge<V, E> edge = adjMatrix.get(index.get(path.get(i))).get(index.get(path.get(i + 1)));
                length += (Integer) edge.weight;
                bandwidth = min(bandwidth, (Integer) edge.bandwidth);
            }
            System.out.printf("%d %d %d\n", path.size(), length, bandwidth);
            for (Vertex<V> vertex : path) {
                System.out.print(vertex.getName() + " ");
            }
            System.out.println();

        } else {
            System.out.println("IMPOSSIBLE");
        }
    }
}