#ifndef DEFORMATION_HH
#define DEFORMAITON_HH

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <GL/glew.h>
#include <GL/GLU.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "ClosestPoint.hh"
#include <gmm.h>
#include <iostream>
#include <set>
#include <float.h>

#define SUBSAMPLE_RADIUS_SCALAR 5
#define SUBSAMPLE_MULTIPIER_FOR_RENDERING 3
#define OPT_ROTATION_WEIGHT 1.0f
#define OPT_REGULARIZATION_WEIGHT 10.0f
#define OPT_CONSTRAINT_WEIGHT 100.0f
#define MAX_ITERATION_STEP 20

#define DEBUG_MODE true

using namespace std;

class Deformation
{

public:
	typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
	typedef OpenMesh::Vec3f Vec3D;
	typedef gmm::dense_matrix<float> gmmMatrix;
	typedef Eigen::VectorXf gmmVector;
	typedef Eigen::SparseMatrix<float> gmmSparseMatrix;
	typedef Eigen::Triplet<float> Trip;
	typedef std::vector<Trip> Triplets;
	
	struct Node
	{
		Vec3D pos, tmpos, t;
		Mesh::VertexHandle vh;
		gmmMatrix R;
		float Ereg;
		bool IsFixed;
		int idx_start;
		set<int> neighbors;
	};

protected:
	typedef vector<int> Influences;
	typedef vector<float> Weights;
	// TODO: Fix Bugs
	typedef OpenMesh::VPropHandleT<Influences> Influencing_Nodes;
	typedef OpenMesh::VPropHandleT<Weights> Vertex_Weights;
	
public:
	Deformation();
	~Deformation();
	bool init(std::string _filename);
	void GenerateGraph();
	Vec3D TransformVerteX(Mesh::VertexHandle v);

	// Optimization functions. Should only be called after handles and fixed points are specified
	// Compute energy function
	float ComputeE();
	// Compute (dim x 1) vector f, where E = fTf
	void Computef();
	// Compute (dim x 12N) Jacobian of f;
	void ComputeJ();
	// Solve for update vector with a Sparse Cholesky Solver
	void SolveLP();
	// Deform mesh vertices based on deformed graph
	bool DeformMesh();
	bool ConvergenceDetection(float fk, float fkk);

	// Support functions for OpenGL
	std::vector<float> get_mesh_buffer_data();
	std::vector<float> get_embedded_graph_buffer_data();
	int getBufferOffsetFromEGBuffer(float x, float y, float z, float min_dist);
	void setRenderedNodePositionUpdate(float x, float y, float z, int idx);
	// For testing only
	void ManuallySetupHandles(int index, Vec3D pos);


protected:
	Mesh mesh_;
	vector<Node> nodes_;
	vector<int> nodes_subsampled_;
	set<int> handle_list_;
	Influencing_Nodes iNodes;
	Vertex_Weights weights;
	

private:
	float get_average_vertex_distance();
	void prepareHandles();

	Vec3D gmmCustomMult(gmmMatrix R, Vec3D v);
	Vec3D gmmCustomMult(Eigen::Matrix3f R, Vec3D v);
	float averageVertexDistance_;
	int id_counter_;
	std::map<Mesh::VertexHandle, Vec3D> Handles_;
	Eigen::SimplicialLLT<gmmSparseMatrix> solver_;
	gmmSparseMatrix Jacobian_;
	gmmVector v_, result_;
	int dimension_of_r, estimation_of_entries;
	bool IfPerformedSymbolicComposition;
	ClosestPoint cp_interactive;
	
};

#endif