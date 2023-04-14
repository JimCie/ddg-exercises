// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        geometry->vertexIndices[v] = idx;
        ++idx;
    }

    idx = 0;
    for (Edge e : mesh->edges()) {
        geometry->edgeIndices[e] = idx;
        ++idx;

    }

    idx = 0;
    for (Face f : mesh->faces()) {
        geometry->faceIndices[f] = idx;
        ++idx;
    }
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    // TODO : DONE
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());

    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh->nEdges()*2);

    size_t rowIndex;
    size_t colIndex;
    for(Edge e : mesh->edges()){
        rowIndex = e.getIndex();
        colIndex = e.firstVertex().getIndex();
        tripletList.push_back(T(rowIndex,colIndex,1));
        colIndex = e.secondVertex().getIndex();
        tripletList.push_back(T(rowIndex,colIndex,1));
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat; // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO : DONE
    SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh->nFaces()*3);
    size_t rowIndex;
    size_t colIndex;
    for(Face f : mesh->faces()){
        rowIndex = f.getIndex();
        for(Edge e : f.adjacentEdges()){
            colIndex = e.getIndex();
            tripletList.push_back(T(rowIndex, colIndex, 1));
        }
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nVertices());
    for(auto index : subset.vertices){
        vec[index] = 1;
    }

    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nEdges());
    for(auto index : subset.edges){
        vec[index] = 1;
    }

    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO : DONE
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nFaces());
    for(auto index : subset.faces){
        vec[index] = 1;
    }

    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO : DONE
    MeshSubset other{};
    /*The result of the multiplication of the adjacency matrices A1*AO
     * is a matrix contains the contribution of the vertex to a face counted twice (one for each edge).
     * */
    SparseMatrix<size_t> B = A1 * A0;

    /*We use the adjacency matrices*/
    for(auto indexVertex : subset.vertices){
        /*go through the column at indexVertex (eigen default storage is column-major*/
        for(SparseMatrix<size_t>::InnerIterator it(A0,indexVertex); it; ++it){
            other.addEdge(it.row());
        }
        for(SparseMatrix<size_t>::InnerIterator it(B,indexVertex); it; ++it){
            other.addFace(it.row());
        }

    }
    for(auto indexEdge : subset.edges){
        for(SparseMatrix<size_t>::InnerIterator it(A1, indexEdge); it; ++it){
            other.addFace(it.row());
        }
    }

    other.addSubset(subset);
    return other; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    auto closure(subset);

    for(auto e : subset.edges){
        for(auto v : mesh->edge(e).adjacentVertices()){
            closure.vertices.insert(v.getIndex());
        }
    }

    for(auto f : subset.faces){
        for(auto v : mesh->face(f).adjacentVertices()){
            closure.vertices.insert(v.getIndex());
        }
        for(auto e : mesh->face(f).adjacentEdges()){
            closure.edges.insert(e.getIndex());
        }
    }

    return closure;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    MeshSubset closureThenStar = star(closure(subset));
    MeshSubset starThenClosure = closure(star(subset));
    starThenClosure.deleteSubset(closureThenStar);
    return starThenClosure;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    MeshSubset cl = closure(subset);
    return cl.equals(subset); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    MeshSubset intermediateMesh{};
    int degree = -1;
    if(!isComplex(subset))
        return -1;
    if(!subset.vertices.empty())
        degree = 0; //this will change if there are simplices other than vertices
    else
        return -1;

    if(!subset.edges.empty()){
        degree = 1;
        //check if all vertices are within an edge
        //this will ensure from previous point that they are within a face if any
        for(size_t indexEdge : subset.edges){
            for(Vertex v : mesh->edge(indexEdge).adjacentVertices()){
                intermediateMesh.addVertex(v.getIndex());
            }
        }
    }else{
        //it only has vertices thus it is pure of degree 0
        return degree;
    }

    if(!subset.faces.empty())
    {
        degree = 2;
        //check if all edges are within a face
        //populate a set with the edges of the faces and compare later
        for(size_t indexFace : subset.faces){
            for(Edge e : mesh->face(indexFace).adjacentEdges()){
                intermediateMesh.addEdge(e.getIndex());
            }
        }
    }else{
        if(intermediateMesh.vertices == subset.vertices)
            return degree;
    }

    if(intermediateMesh.edges == subset.edges and intermediateMesh.vertices == subset.vertices)
        return degree;
    return -1; // if all previous tests have failed, which should not happen
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO : DONE
    MeshSubset preBoundary;
    size_t contained_in = 0;
    /*if a vertex occurs in only one edge, it is a boundary*/
    /*we check if there are many of the adjacent edges in the subset*/
    for(auto indexVertex : subset.vertices){
        for(Edge e : mesh->vertex(indexVertex).adjacentEdges()){
            if(subset.edges.count(e.getIndex())){
                contained_in++;
            }
        }
        if(contained_in == 1){
            preBoundary.addVertex(indexVertex);
        }
        contained_in = 0;
    }

    /*if an edge occurs in only one face, it is a boundary*/
    for(auto indexEdge : subset.edges){
        for(Face f : mesh->edge(indexEdge).adjacentFaces()){
            if(subset.faces.count(f.getIndex())){
                contained_in++;
            }
        }
        if(contained_in == 1){
            preBoundary.addEdge(indexEdge);
        }
        contained_in = 0;
    }

    return closure(preBoundary); // placeholder
}