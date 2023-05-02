// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    SparseMatrix<double> star0(mesh.nVertices(), mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nVertices());
    /* //just for fun
    for(size_t index = 0; index<mesh.nVertices(); index++){
        Vertex v = mesh.vertex(index);
        double area = VertexPositionGeometry::barycentricDualArea(v); // no need to divide by the measure of the vertex, we take it to 1
        tripletList.push_back(T(index, index, area));
    }*/

    for(Vertex v : mesh.vertices()){ //other possible traversal
        double area = VertexPositionGeometry::barycentricDualArea(v); // no need to divide by the measure of the vertex, we take it to 1
        tripletList.push_back(T(v.getIndex(), v.getIndex(), area));
    }

    star0.setFromTriplets(tripletList.begin(),  tripletList.end());
    return star0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;

    for(auto e : mesh.edges()){

        double cot1 = VertexPositionGeometry::cotan(e.halfedge());
        double cot2 = VertexPositionGeometry::cotan(e.halfedge().twin());

        double value = (cot1 + cot2) / 2;
        triplets.push_back(T(e.getIndex(), e.getIndex(), value ));
    }

    Eigen::SparseMatrix<double> sparseMatrix(mesh.nEdges(), mesh.nEdges());
    sparseMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return sparseMatrix;
}
/*SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    SparseMatrix<double> star1(mesh.nEdges(), mesh.nEdges());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nEdges());

    for(Edge e : mesh.edges()){
        double value;
        if(e.isBoundary()){
            value = VertexPositionGeometry::cotan(e.halfedge());
        }else{
            value = VertexPositionGeometry::cotan(e.halfedge()) + VertexPositionGeometry::cotan(e.halfedge().twin());
        }

        value *= 0.5;
        value/= edgeLength(e);
        tripletList.push_back(T(e.getIndex(), e.getIndex(), value));
    }

    star1.setFromTriplets(tripletList.begin(), tripletList.end());
    return star1;
}*/

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    SparseMatrix<double> star2(mesh.nFaces(), mesh.nFaces());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nFaces());

    for(auto f: mesh.faces()){ // f Face
        tripletList.push_back(T(f.getIndex(), f.getIndex(), double(1./faceArea(f))));
    }

    star2.setFromTriplets(tripletList.begin(), tripletList.end());
    return star2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    SparseMatrix<double> d0(mesh.nEdges(), mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nEdges()*2);
    for(auto e : mesh.edges()){
        auto tail = e.firstVertex().getIndex();
        auto tip = e.secondVertex().getIndex();
        tripletList.push_back(T(e.getIndex(), tip, 1));
        tripletList.push_back(T(e.getIndex() , tail, -1));
    }
    d0.setFromTriplets(tripletList.begin(), tripletList.end());
    return d0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    SparseMatrix<double> d1(mesh.nFaces(), mesh.nEdges());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(3*mesh.nFaces());
    for(auto f : mesh.faces()){
        for(auto he : f.adjacentHalfedges()){
            double orientation = he.orientation() ? 1 : -1;
            tripletList.push_back(T(f.getIndex(), he.edge().getIndex(), orientation));
        }
    }
    d1.setFromTriplets(tripletList.begin(), tripletList.end());
    return d1;
}

} // namespace surface
} // namespace geometrycentral