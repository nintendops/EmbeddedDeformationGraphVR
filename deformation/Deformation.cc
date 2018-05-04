#include "Deformation.hh"

// Constructor
Deformation::Deformation()
{
}

// Destructor
Deformation::~Deformation()
{

}

bool Deformation::init(std::string _filename)
{
	// initialization
	std::cout << "Reading:  " << _filename << std::endl;

	nodes_.clear();
	handle_list_.clear();
	dimension_of_r = 0;
	estimation_of_entries = 0;
	mesh_.request_vertex_status();
	mesh_.request_vertex_normals();
	mesh_.request_face_normals();
	mesh_.request_edge_status();
	mesh_.request_face_status();

	if (OpenMesh::IO::read_mesh(mesh_, _filename)) 
	{
		std::cout << "#vertices: " << mesh_.n_vertices() << std::endl;
		mesh_.update_normals();
		mesh_.add_property(iNodes);
		mesh_.add_property(weights);
		
		if (!mesh_.has_vertex_normals())
		{
			std::cout << "No vertex normal read!"<< std::endl;
			return false;
		}

		// scaling up the mesh a bit
		float scale = 5.0f;
		Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
		for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
		{
			Vec3D p = scale * mesh_.point(*v_it);
			mesh_.set_point(*v_it, p);
		}
		averageVertexDistance_ = get_average_vertex_distance();

		return true;
	}
	else
		return false;
}

void Deformation::GenerateGraph()
{
	Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
	Mesh::VertexHandle v1,v2;
	vector<Node>::iterator vn_begin, vn_end;
	vector<Vec3D> npos;
	//ClosestPoint cp;
	Vec3D p1, p2, pv;
	int node_count = 0;

	// Generate embedded nodes by uniform sampling
	set<OpenMesh::VertexHandle> candidates;
	float subsampleRadius = SUBSAMPLE_RADIUS_SCALAR * averageVertexDistance_;
	for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)		candidates.insert(*v_it);
	while (!candidates.empty())
	{
		v1 = *candidates.begin();
		p1 = mesh_.point(v1);
		// add node to embedded graph
		Node n;
		n.pos = p1;
		n.tmpos = p1;
		//n.normal = mesh_.normal(v1);
		n.IsFixed = false;
		n.vh = v1;
		//Hackhack. Should change this if you have fixed handle
		n.idx_start = node_count * 12;
		nodes_.push_back(n);
		node_count++;

		// remove nodes that are close
		set<OpenMesh::VertexHandle>::iterator s_it = candidates.begin();
		for (; s_it != candidates.end();)
		{
			v2 = *s_it;
			p2 = mesh_.point(v2);
			if (sqrt((p1 - p2) | (p1 - p2)) < subsampleRadius)
			{
				candidates.erase(s_it++);
				if (candidates.size() == 0)
					s_it = candidates.end();
			}
			else
				++s_it;
		}
	}

	std::cout << "sampled #" << nodes_.size() << " EG nodes" << std::endl;

	// set up ann structure
	vn_begin = nodes_.begin();
	vn_end = nodes_.end();
	for(; vn_begin != vn_end; ++vn_begin)		npos.push_back(vn_begin->pos);
	cp_interactive.init(npos);

	std::cout << "finding influence nodes and neighbors. Also computing weights..." << std::endl;

	for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
	{
		// finding k nearest nodes 
		pv = mesh_.point(*v_it);
		vector<int> bestIdx = cp_interactive.getClosestPoints(pv, K_NEAREST_NODE + 1);
		for(int j = 0; j < K_NEAREST_NODE; j++) 	mesh_.property(iNodes, *v_it).push_back(bestIdx[j]);

		// constructing graph edges
		for (int j = 0; j < K_NEAREST_NODE; j++)
		{
			for (int jj = 0; jj < K_NEAREST_NODE; jj++)
			{
				if (jj != j)
				{
					// every edge gives 3 regularization term (x,y,z) with 12 non zero entries
					dimension_of_r += 3;
					estimation_of_entries += 12;
					nodes_[bestIdx[j]].neighbors.insert(bestIdx[jj]);
				}
			}
		}

		// computnig weights
		Vec3D pmax = npos[bestIdx[K_NEAREST_NODE]];		
		float dmax = sqrt((pv - pmax) | (pv - pmax));
		float weight_sum = 0.0f;
		float weight[K_NEAREST_NODE];
		for (int j = 0; j < K_NEAREST_NODE; j++)
		{
			Vec3D g = npos[bestIdx[j]];
			float pg_dist = sqrt((pv - g) | (pv - g));
			weight[j] = pow((1.0f - pg_dist / dmax), 2);
			weight_sum += weight[j];
		}
		float scale = 1 / weight_sum;
		for (int j = 0; j < K_NEAREST_NODE; j++)
			mesh_.property(weights, *v_it).push_back(weight[j] * scale);
	}

	std::cout << "Initializing variables..." << std::endl;

	// Initialize R,t per node
	// R is initialized to be slightly different from identity I. t initialized to be zero
	gmmMatrix init_r(3, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			init_r(i, j) = 0.1f;
	init_r(0, 0) = 0.9f;
	init_r(1, 1) = 0.9f;
	init_r(2, 2) = 0.9f;
	Vec3D init_t(0,0,0);

	for(vn_begin = nodes_.begin(); vn_begin != vn_end; ++vn_begin)
	{
		// every R gives 6 rotation terms with 27 non_zero entries
		dimension_of_r += 6;
		estimation_of_entries += 27;
		vn_begin->R = init_r;
		vn_begin->t = init_t;
	}

	// further subsample graph nodes for rendering
	set<int> nodes_candidates;
	int i1, i2;
	Vec3D n1, n2;
	subsampleRadius *= SUBSAMPLE_MULTIPIER_FOR_RENDERING;
	for (int i = 0; i < nodes_.size(); i++) nodes_candidates.insert(i);

	while (!nodes_candidates.empty())
	{
		i1 = *nodes_candidates.begin();
		n1 = nodes_[i1].pos;
		nodes_subsampled_.push_back(i1);

		// remove nodes that are close
		set<int>::iterator s_it = nodes_candidates.begin();
		for (; s_it != nodes_candidates.end();)
		{
			i2 = *s_it;
			n2 = nodes_[i2].pos;

			if (sqrt((n1 - n2) | (n1 - n2)) < subsampleRadius)
			{
				nodes_candidates.erase(s_it++);
				if (nodes_candidates.size() == 0)
					s_it = nodes_candidates.end();
			}
			else
				++s_it;
		}
	}

	std::cout << "Further subsampled #" << nodes_subsampled_.size() << " nodes for rendering" << std::endl;
	cp_interactive.release();
}


float Deformation::get_average_vertex_distance()
{
	float accDist = 0;
	int accCount = 0;

	Mesh::ConstHalfedgeIter he_it = mesh_.halfedges_begin();
	for (; he_it != mesh_.halfedges_end(); ++he_it)
	{
		Vec3D p = mesh_.point(mesh_.from_vertex_handle(*he_it));
		Vec3D q = mesh_.point(mesh_.to_vertex_handle(*he_it));
		float edgeLength = sqrt((p - q) | (p - q));
		accDist += edgeLength;
		accCount++;
	}

	if (accCount > 0) return accDist / float(accCount);
	else return 0.0f;
}

Deformation::Vec3D Deformation::TransformVerteX(Mesh::VertexHandle v)
{
	Node j;
	Vec3D v1;
	float w;
	Vec3D vi = mesh_.point(v);
	Vec3D v2(0, 0, 0); 

	for (int i = 0; i < K_NEAREST_NODE; i++)
	{
		j = nodes_[mesh_.property(iNodes, v)[i]];
		w = mesh_.property(weights, v)[i];
		v1 = gmmCustomMult(j.R, vi - j.pos);
		v1 = w * (v1 + j.pos + j.t);
		v2 += v1;
	}

	return v2;
}

float Deformation::ComputeE()
{
	std::vector<Deformation::Node>::iterator n_it = nodes_.begin();
	std::set<int>::iterator nn_it;
	gmmMatrix R;
	Vec3D tk, gk, gj, tj, v1;

	float en = 0.0f;
	float en_rot = 0.0f;
	float en_reg = 0.0f;
	float en_con = 0.0f;

	for (; n_it != nodes_.end(); ++n_it)
	{
		if (n_it->IsFixed) continue;

		R = n_it->R;
		tk = n_it->t;
		gk = n_it->pos;

		// Erot
		Vec3D c1(R(0, 0), R(1, 0), R(2, 0));
		Vec3D c2(R(0, 1), R(1, 1), R(2, 1));
		Vec3D c3(R(0, 2), R(1, 2), R(2, 2));
		//en += OPT_ROTATION_WEIGHT * (pow(dot(c1, c2), 2) + pow(dot(c1, c3), 2) + pow(dot(c2, c3), 2)
		//	+ pow(dot(c1, c1) - 1, 2) + pow(dot(c2, c2) - 1, 2) + pow(dot(c3, c3) - 1, 2));

		en_rot += OPT_ROTATION_WEIGHT * (pow(dot(c1, c2), 2) + pow(dot(c1, c3), 2) + pow(dot(c2, c3), 2)
			+ pow(dot(c1, c1) - 1, 2) + pow(dot(c2, c2) - 1, 2) + pow(dot(c3, c3) - 1, 2));

		// Ereg
		for (nn_it = n_it->neighbors.begin(); nn_it != n_it->neighbors.end(); ++nn_it)
		{
			gj = nodes_[*nn_it].pos;
			tj = nodes_[*nn_it].t;
			v1 = gmmCustomMult(R, (gk - gj));
			v1 = v1 + gj + tj - (gk + tk);
			//en += OPT_REGULARIZATION_WEIGHT * (v1 | v1);
			en_reg += OPT_REGULARIZATION_WEIGHT * (v1 | v1);
		}
	}

	// Econ
	std::map<Mesh::VertexHandle, Vec3D>::iterator vh_it = Handles_.begin();
	Mesh::VertexHandle vh;
	Vec3D ql, qd;

	for (; vh_it != Handles_.end(); ++vh_it)
	{
		vh = vh_it->first;
		ql = vh_it->second;
		qd = TransformVerteX(vh) - ql;
		//en += OPT_CONSTRAINT_WEIGHT * (qd | qd);
		en_con += OPT_CONSTRAINT_WEIGHT * (qd | qd);
	}
	en = en_rot + en_reg + en_con;
	if (DEBUG_MODE)
	{
		std::cout << "Computed E rot at this iteration: " << en_rot << std::endl;
		std::cout << "Computed E reg at this iteration: " << en_reg << std::endl;
		std::cout << "Computed E con at this iteration: " << en_con << std::endl;
		std::cout << "Computed E at this iteration: " << en << std::endl;
	}
	return en;
}

void Deformation::Computef()
{
	//std::cout << "Computing f" << std::endl;

	v_ = gmmVector(dimension_of_r);
	v_.setZero();
	gmmMatrix R;
	std::vector<Deformation::Node>::iterator n_it = nodes_.begin();
	std::set<int>::iterator nn_it;
	Vec3D t, tk, gk, gj, v1;
	int idx = 0;

	for (; n_it != nodes_.end(); ++n_it)
	{
		if (n_it->IsFixed) continue;
		R = n_it->R;
		t = n_it->t;
		gj = n_it->pos;
		Vec3D c1(R(0, 0), R(1, 0), R(2, 0));
		Vec3D c2(R(0, 1), R(1, 1), R(2, 1));
		Vec3D c3(R(0, 2), R(1, 2), R(2, 2));

		// rotational terms
		v_[idx] = (sqrt(OPT_ROTATION_WEIGHT)*(dot(c1, c2)));
		v_[idx + 1] = (sqrt(OPT_ROTATION_WEIGHT)*(dot(c1, c3)));
		v_[idx + 2] = (sqrt(OPT_ROTATION_WEIGHT)*(dot(c2, c3)));
		v_[idx + 3] = (sqrt(OPT_ROTATION_WEIGHT)*(dot(c1, c1) - 1));
		v_[idx + 4] = (sqrt(OPT_ROTATION_WEIGHT)*(dot(c2, c2) - 1));
		v_[idx + 5] = (sqrt(OPT_ROTATION_WEIGHT)*(dot(c3, c3) - 1));
		idx += 6;

		// regularization terms
		for (nn_it = n_it->neighbors.begin(); nn_it != n_it->neighbors.end(); ++nn_it)
		{
			gk = nodes_[*nn_it].pos;
			tk = nodes_[*nn_it].t;
			v1 = gmmCustomMult(R, gk - gj);
			v1 = v1 + gj + t - (gk + tk);
			v_[idx] = (sqrt(OPT_REGULARIZATION_WEIGHT)*v1[0]);
			v_[idx + 1] = (sqrt(OPT_REGULARIZATION_WEIGHT)*v1[1]);
			v_[idx + 2] = (sqrt(OPT_REGULARIZATION_WEIGHT)*v1[2]);
			idx += 3;
		}
	}

	// constraint terms
	std::map<Mesh::VertexHandle, Vec3D>::iterator vh_it = Handles_.begin();
	Mesh::VertexHandle vh;
	Vec3D ql, qd;

	for (; vh_it != Handles_.end(); ++vh_it)
	{
		vh = vh_it->first;
		ql = vh_it->second;
		qd = TransformVerteX(vh) - ql;
		v_[idx] = (sqrt(OPT_CONSTRAINT_WEIGHT)*qd[0]);
		v_[idx + 1] = (sqrt(OPT_CONSTRAINT_WEIGHT)*qd[1]);
		v_[idx + 2] = (sqrt(OPT_CONSTRAINT_WEIGHT)*qd[2]);
		idx += 3;
	}

	// Testing only
	float en = 0.0f;
	for (int i = 0; i < v_.size(); i++)
		en += v_[i] * v_[i];
	if(DEBUG_MODE)
		std::cout << "Computed fTf at this iteration: " << en << std::endl;

}

void Deformation::ComputeJ()
{
	//std::cout << "Computing J" << std::endl;
	int n_size = nodes_.size();
	Jacobian_ = Eigen::SparseMatrix<float>(dimension_of_r, 12 * n_size);
	Jacobian_.setZero();
	Triplets tList;
	tList.reserve(estimation_of_entries);
	gmmMatrix R;
	std::vector<Deformation::Node>::iterator n_it = nodes_.begin();
	std::set<int>::iterator nn_it;
	Vec3D t, tk, gk, gj, v1;

	// k1: index counter for row; 
	// k2 : index counter for column (starting index of variables for node i)
	int k1 = 0, k2 = 0;

	float wwRot = sqrt(OPT_ROTATION_WEIGHT);
	float wwReg = sqrt(OPT_REGULARIZATION_WEIGHT);

	for (; n_it != nodes_.end(); ++n_it)
	{
		if (n_it->IsFixed) continue;
		R = n_it->R;
		t = n_it->t;
		gj = n_it->pos;

		// record starting index of the node variables in row vector (probably should just precompute this)
		// n_it->idx_start = k2;

		// Erot 
		// c1c2
		for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1, 3 * i + k2, R(i, 1) * wwRot)); tList.push_back(Trip(k1, 3 * i + k2 + 1, R(i, 0) * wwRot)); }
		// c1c3
		for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 1, 3 * i + k2, R(i, 2) * wwRot)); tList.push_back(Trip(k1 + 1, 3 * i + k2 + 2, R(i, 0) * wwRot)); }
		// c2c3
		for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 2, 3 * i + k2 + 1, R(i, 2) * wwRot)); tList.push_back(Trip(k1 + 2, 3 * i + k2 + 2, R(i, 1) * wwRot)); }
		// cc - 1
		for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 3, 3 * i + k2, 2 * R(i, 0) * wwRot));}
		for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 4, 3 * i + k2 + 1, 2 * R(i, 1) * wwRot));}
		for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 5, 3 * i + k2 + 2, 2 * R(i, 2) * wwRot));}
		k1 += 6;

		// Ereg
		for (nn_it = n_it->neighbors.begin(); nn_it != n_it->neighbors.end(); ++nn_it)
		{
			// tk = gk - gj
			int gk_idx = nodes_[*nn_it].idx_start;
			tk = nodes_[*nn_it].pos - gj;
			for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1, k2 + i, wwReg * tk[i])); }
			for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 1, k2 + 3 + i, wwReg * tk[i])); }
			for (int i = 0; i < 3; i++) { tList.push_back(Trip(k1 + 2, k2 + 6 + i, wwReg * tk[i])); }
			tList.push_back(Trip(k1, k2 + 9, wwReg));
			tList.push_back(Trip(k1 + 1, k2 + 10, wwReg));
			tList.push_back(Trip(k1 + 2, k2 + 11, wwReg));
			if (!nodes_[*nn_it].IsFixed)
			{
				tList.push_back(Trip(k1, gk_idx + 9, -wwReg));
				tList.push_back(Trip(k1 + 1, gk_idx + 10, -wwReg));
				tList.push_back(Trip(k1 + 2, gk_idx + 11, -wwReg));
			}
			k1 += 3;
		}
		k2 += 12;
	}

	//Econ
	std::map<Mesh::VertexHandle, Vec3D>::iterator vh_it = Handles_.begin();
	Mesh::VertexHandle vh;
	Vec3D ql, qd;
	Node j;
	float w;
	int idx;

	for (; vh_it != Handles_.end(); ++vh_it)
	{
		vh = vh_it->first;
		ql = vh_it->second;

		for (int i = 0; i < K_NEAREST_NODE; i++)
		{
			j = nodes_[mesh_.property(iNodes, vh)[i]];
			if (j.IsFixed) continue;

			w = mesh_.property(weights, vh)[i] * sqrt(OPT_CONSTRAINT_WEIGHT);
			qd = mesh_.point(vh) - j.pos;
			idx = j.idx_start;

			// each handle constraint adds three rows in J
			for (int ii = 0; ii < 3; ii++) tList.push_back(Trip(k1, idx + ii, w*qd[ii]));
			for (int ii = 0; ii < 3; ii++) tList.push_back(Trip(k1 + 1, idx + ii + 3, w*qd[ii]));
			for (int ii = 0; ii < 3; ii++) tList.push_back(Trip(k1 + 2, idx + ii + 6, w*qd[ii]));
			tList.push_back(Trip(k1, idx + 9, w));
			tList.push_back(Trip(k1 + 1, idx + 10, w));
			tList.push_back(Trip(k1 + 2, idx + 11, w));
		}
		k1 += 3;
	}
	Jacobian_.setFromTriplets(tList.begin(), tList.end());
}

// Solve for update vector with a Sparse Cholesky Solver
void Deformation::SolveLP()
{
	gmmSparseMatrix A = Jacobian_.transpose() * Jacobian_;
	gmmVector b = -1.0f * Jacobian_.transpose() * v_;
	
	if (!IfPerformedSymbolicComposition) 
	{
		solver_.analyzePattern(A); 
		IfPerformedSymbolicComposition = true;
	}
	solver_.factorize(A);
	result_ = solver_.solve(b);

	// update variables
	std::vector<Deformation::Node>::iterator n_it = nodes_.begin();
	gmmMatrix Rk;
	Vec3D tk;

	for (; n_it != nodes_.end(); ++n_it)
	{
		if (n_it->IsFixed)
			continue;

		int idx = n_it->idx_start;
		Rk = n_it->R;
		tk = n_it->t;
		Rk(0, 0) += result_(idx);
		Rk(0, 1) += result_(idx + 1);
		Rk(0, 2) += result_(idx + 2);
		Rk(1, 0) += result_(idx + 3);
		Rk(1, 1) += result_(idx + 4);
		Rk(1, 2) += result_(idx + 5);
		Rk(2, 0) += result_(idx + 6);
		Rk(2, 1) += result_(idx + 7);
		Rk(2, 2) += result_(idx + 8);
		tk[0] += result_(idx + 9);
		tk[1] += result_(idx + 10);
		tk[2] += result_(idx + 11);
		n_it->R = Rk;
		n_it->t = tk;
	}
	
}

void Deformation::prepareHandles()
{
	set<int>::iterator s_it = handle_list_.begin();
	Mesh::VertexHandle vh;

	for (; s_it != handle_list_.end(); s_it++)
	{
		dimension_of_r += 3;
		Handles_.insert(std::pair<Mesh::VertexHandle, Vec3D>(nodes_[*s_it].vh, nodes_[*s_it].tmpos));
	}
}

bool Deformation::DeformMesh()
{
	if (handle_list_.empty())
		return false;

	prepareHandles();
	IfPerformedSymbolicComposition = false;
	if (DEBUG_MODE) 
	{
		std::cout << "Start deforming mesh. Number of variables: " << nodes_.size() * 12
			<< "; Number of row entries: " << dimension_of_r << std::endl;
	}
	Mesh::VertexIter v_it = mesh_.vertices_begin();
	Vec3D v0, v1, n0, gj;
	Node nj;
	float fk, fkk, wj;

	int iteration_counter = 0;
	fk = ComputeE();
	do 
	{
		fkk = fk;
		Computef();
		ComputeJ();
		SolveLP();
		fk = ComputeE();
		iteration_counter++;
	} while (!ConvergenceDetection(fk, fkk) && iteration_counter < MAX_ITERATION_STEP);

	std::cout << "# of iteration taken to optimize: " << iteration_counter << std::endl;

	std::cout << "Deforming mesh vertices..." << std::endl;

	// Perform "Skeleton" Subspace Deformation to mesh vertices	
	int dc = 0;

	for (; v_it != mesh_.vertices_end(); ++v_it)
	{
		v0 = mesh_.point(*v_it);
		n0 = mesh_.normal(*v_it);
		Vec3D vk(0, 0, 0);
		Vec3D nk(0, 0, 0);
		for (int j = 0; j < K_NEAREST_NODE; j++)
		{
			wj = mesh_.property(weights, *v_it)[j];
			nj = nodes_[mesh_.property(iNodes, *v_it)[j]];

			vk += wj * (gmmCustomMult(nj.R, v0 - nj.pos) + nj.pos + nj.t);

			// stupid gmm not let me do inverse...HACKHACK
			Eigen::Matrix3f RR;
			RR << nj.R(0, 0), nj.R(0, 1), nj.R(0, 2),
				  nj.R(1, 0), nj.R(1, 1), nj.R(1, 2),
				  nj.R(2, 0), nj.R(2, 1), nj.R(2, 2);
			nk += wj * gmmCustomMult(RR.inverse().transpose(), n0);
		}
		if (DEBUG_MODE && dc < 4)
		{
			std::cout << "Vertex position before deformation:" << v0 << std::endl;
			std::cout << "Vertex position after deformation:" << vk << std::endl;
			dc++;
		}

		mesh_.set_point(*v_it, vk);
		mesh_.set_normal(*v_it, nk);
	}
	mesh_.update_normals();

	// Move the position of the embedded graph nodes too
	std::vector<Deformation::Node>::iterator n_it = nodes_.begin();
	for (; n_it != nodes_.end(); ++n_it)
	{
		n_it->pos += n_it->t;
		n_it->tmpos = n_it->pos;
	}

	// Reinitializing stuffs
	Handles_.clear();
	handle_list_.clear();
	return true;
}


bool Deformation::ConvergenceDetection(float fk, float fkk)
{
	float epsilon = 1e-06F;
	bool objCond = abs(fk - fkk) < epsilon * (1 + fk);
	// TODO: gradient condition

	float maxUpdate = -1.0f;
	for (int i = 0; i < result_.size(); i++)
	{
		if (abs(result_(i)) > maxUpdate) maxUpdate = abs(result_(i));
	}
	bool updateCond = (maxUpdate < sqrt(epsilon) * (1 + maxUpdate));
	return objCond && updateCond;
}

Deformation::Vec3D Deformation::gmmCustomMult(gmmMatrix R, Vec3D v)
{
	Vec3D Rv;
	Rv[0] = R(0, 0) * v[0] + R(0, 1) * v[1] + R(0, 2) * v[2];
	Rv[1] = R(1, 0) * v[0] + R(1, 1) * v[1] + R(1, 2) * v[2];
	Rv[2] = R(2, 0) * v[0] + R(2, 1) * v[1] + R(2, 2) * v[2];
	return Rv;
}

Deformation::Vec3D Deformation::gmmCustomMult(Eigen::Matrix3f R, Vec3D v)
{
	Vec3D Rv;
	Rv[0] = R(0, 0) * v[0] + R(0, 1) * v[1] + R(0, 2) * v[2];
	Rv[1] = R(1, 0) * v[0] + R(1, 1) * v[1] + R(1, 2) * v[2];
	Rv[2] = R(2, 0) * v[0] + R(2, 1) * v[1] + R(2, 2) * v[2];
	return Rv;
}

//OpenGL support functions
std::vector<float> Deformation::get_mesh_buffer_data()
{
	std::vector<float> buffer;
	Mesh::FaceIter f_it = mesh_.faces_begin();
	Mesh::FVIter fv_it;
	Vec3D v,normal;
	for (; f_it != mesh_.faces_end(); ++f_it)
	{
		// should make sure faces are triangle
		int i = 0;
		for (fv_it = mesh_.fv_begin(*f_it); fv_it.is_valid(); ++fv_it)
		{
			i++;
			if (i > 3)
			{
				std:cerr << "Mesh has non-triangle face!";
				break;
			}

			v = mesh_.point(*fv_it);
			normal = mesh_.normal(*fv_it);

			//pos
			buffer.push_back(v[0]);
			buffer.push_back(v[1]);
			buffer.push_back(v[2]);
			//normal
			buffer.push_back(normal[0]);
			buffer.push_back(normal[1]);
			buffer.push_back(normal[2]);
		}
	}

	return buffer;
}

// does not include information about connectivity
std::vector<float> Deformation::get_embedded_graph_buffer_data()
{
	std::vector<float> buffer;
	std::vector<int>::iterator n_it = nodes_subsampled_.begin();
	for (; n_it != nodes_subsampled_.end(); ++n_it)
	{
		buffer.push_back(nodes_[*n_it].pos[0]);
		buffer.push_back(nodes_[*n_it].pos[1]);
		buffer.push_back(nodes_[*n_it].pos[2]);
	}
	return buffer;
}

// return -1 if closest node is too far
int Deformation::getBufferOffsetFromEGBuffer(float x, float y, float z, float min_dist)
{
	Vec3D cpos(x, y, z);

	if (!cp_interactive.IfInitialized)
	{
		vector<Vec3D> npos;
		std::vector<int>::iterator n_it = nodes_subsampled_.begin();
		for (; n_it != nodes_subsampled_.end(); ++n_it) npos.push_back(nodes_[*n_it].tmpos);
		cp_interactive.init(npos);
	}

	int bestIdx = cp_interactive.getClosestPoint(cpos);
	Vec3D displacement = cpos - nodes_[nodes_subsampled_[bestIdx]].tmpos;
	if (sqrt(displacement | displacement) <= min_dist * averageVertexDistance_ * 5)
	{
		Vec3D np(x, y, z);
		nodes_[nodes_subsampled_[bestIdx]].tmpos = np;
		cp_interactive.release();
		handle_list_.insert(nodes_subsampled_[bestIdx]);
		return bestIdx;
	}
	else return -1;
}

void Deformation::setRenderedNodePositionUpdate(float x, float y, float z, int idx)
{
	if (idx >= 0 && idx < nodes_subsampled_.size())
	{
		nodes_[nodes_subsampled_[idx]].tmpos = Vec3D(x, y, z);
		handle_list_.insert(nodes_subsampled_[idx]);
	}
	else
		std::cerr << "Wrong index passed into setRenderedNodePositionUpdate()" << std::endl;
}


void Deformation::ManuallySetupHandles(int index, Vec3D pos)
{
	int true_idx = nodes_subsampled_[index];
	nodes_[true_idx].tmpos = pos;
	handle_list_.insert(true_idx);	
}