// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

        size_t v1 = he.tipVertex().getIndex();
        size_t v2 = he.tailVertex().getIndex();

        size_t v3 = he.next().tipVertex().getIndex();

        assert(v1!=v2);
        assert(v1!=v3);
        assert(v2!=v3);

        auto u = vertexPositions[v1] - vertexPositions[v3];
        auto v = vertexPositions[v2] - vertexPositions[v3];

        double denominator = cross(u,v).norm();
        assert(0!=denominator);

        return dot(u,v) / denominator;
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    double area = 0.;
    for(auto f : v.adjacentFaces()){
        area += faceArea(f);
    }
    return area/3.; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    auto v1 = c.vertex().getIndex();
    auto v2 = c.halfedge().tipVertex().getIndex();
    auto v3 = c.halfedge().next().tipVertex().getIndex();

    const auto u = vertexPositions[v2] - vertexPositions[v1];
    const auto v = vertexPositions[v3] - vertexPositions[v1];

    double denominator = norm(u)*norm(v);
    assert(0!= denominator);

    return clamp(acos(dot(u,v) / denominator), 0.0, PI);
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
    auto n1 = faceNormals[he.face().getIndex()];//must add geometry->requireFaceNormals();
    auto n2 = faceNormals[he.twin().face().getIndex()];

    auto E = vertexPositions[he.tipVertex().getIndex()] - vertexPositions[he.tailVertex().getIndex()];
    auto e = E.normalize();

    auto c = cross(n1,n2);
    auto d = dot(n1,n2);

    return atan2(dot(e, c), dot(n1,n2));
/*    Halfedge heNext = he.next(); //CA face 1 : ABC
    Halfedge heTwinNext = he.twin().next(); //BD face 2 : BCD

    const auto BC = vertexPositions[he.tipVertex().getIndex()] - vertexPositions[he.tailVertex().getIndex()];
    const auto CB = -1 * BC;
    const auto CA = vertexPositions[heNext.tipVertex().getIndex()] - vertexPositions[heNext.tailVertex().getIndex()];
    const auto BD = vertexPositions[heTwinNext.tipVertex().getIndex()] - vertexPositions[heTwinNext.tailVertex().getIndex()];

    const auto normal1 = cross(CA, CB).normalize();
    const auto normal2 = cross(BD, BC).normalize();

    const auto edge = BC.normalize();

    return atan2(dot(edge,cross(normal1, normal2)), dot(normal1, normal2)); // placeholder*/
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    auto normal = Vector3::zero();
    for(const auto & f : v.adjacentFaces()){
        normal += faceNormals[f.getIndex()];
    }
    return normal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    auto normal = Vector3::zero();
    for(auto c : v.adjacentCorners()){
        normal += angle(c) * faceNormals[c.face().getIndex()];
    }
    return normal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    auto normal = Vector3::zero();
    auto he1 = v.halfedge();
    auto he2 = he1.twin().next();

    //Vertex is A, we compute the angle for triangles ABC
    while(1){
        Vector3 AB = vertexPositions[he1.tipVertex().getIndex()] - vertexPositions[he1.tailVertex().getIndex()];
        Vector3 AC = vertexPositions[he2.tipVertex().getIndex()] - vertexPositions[he2.tailVertex().getIndex()];
        normal += cross(AC, AB) / (norm2(AB) * norm2(AC));
        he1 = he2;
        he2 = he1.twin().next();
        if(he1 == v.halfedge())
            break;
    }
    return normal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    auto normal = Vector3::zero();
    for(auto f : v.adjacentFaces()){
        normal += faceNormals[f.getIndex()]*faceAreas[f.getIndex()]; // requireFaceAreas();
        //normal += faceNormal(f) * faceArea(f);
    }
    return normal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    auto normal = Vector3::zero();

    for(auto he : v.outgoingHalfedges()){
        const auto posVec = (vertexPositions[he.tipVertex().getIndex()] - vertexPositions[he.tailVertex().getIndex()]).normalize();
        normal += dihedralAngle(he) * posVec;
    }
    return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    auto normal = Vector3::zero();
    for(auto he : v.outgoingHalfedges()){
        const auto posVec = vertexPositions[he.tipVertex().getIndex()] - vertexPositions[he.tailVertex().getIndex()];
        normal += (cotan(he) + cotan(he.twin())) * posVec;
    }
    return normal.normalize();
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    double totalInnerAngle = 0.0;

    for(const auto& c : v.adjacentCorners()){
        totalInnerAngle += angle(c);
    }

    return 2*PI - totalInnerAngle;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    double total = 0.0;

    for(const auto& v : mesh.vertices()){
        total += angleDefect(v);
    }
    return total;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    double meanCurvature = 0.;
    for(auto he : v.outgoingHalfedges()){
        meanCurvature += dihedralAngle(he) * norm(vertexPositions[he.tipVertex().getIndex()]- vertexPositions[he.tailVertex().getIndex()]);
    }
    return meanCurvature * 0.5; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    double area = 0.;
    for(auto he : v.outgoingHalfedges()){
        area += norm2(vertexPositions[he.tipVertex().getIndex()]- vertexPositions[he.tailVertex().getIndex()]) * cotan(he);
        he = he.next().next();
        area += norm2(vertexPositions[he.tipVertex().getIndex()]- vertexPositions[he.tailVertex().getIndex()]) * cotan(he);
    }
    return area * 0.125; // placeholder
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    return std::make_pair(0, 0); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral