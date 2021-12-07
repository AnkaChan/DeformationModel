#include <direct.h>
#include <set>

#include <iostream>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/rotation.h"

#include <MeshFrame/core/Mesh/BaseMesh.h>
#include <MeshFrame/core/Mesh/Vertex.h>
#include <MeshFrame/core/Mesh/HalfEdge.h>
#include <MeshFrame/core/Mesh/Edge.h>
#include <MeshFrame/core/Mesh/Face.h>
#include <MeshFrame/core/Mesh/Types.h>
#include <MeshFrame/core/Mesh/Iterators2.h>

#include <ctime>
#include <stdio.h>
#include "AC/Parser.h"
#include "AC/IO_AC.h"

using namespace MeshLib;

using std::cout;
using std::endl;

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

typedef CBaseMesh<CVertex, CEdge, CFace, CHalfEdge> M;
typedef CIterators<M> It;

struct Config
{
	std::string inExtName = "ply";
	int maxIteration = 100;
	double functionTolerance = 1e-8;

	int startFrameId = 0;
	int endFrameId = -1;
};

struct L2Regularizer
{
	L2Regularizer(double px, double py) :
		px_(px), py_(py)
	{}
	template <typename T> static void initailize(Eigen::Matrix<T, 3, 1> & vec, T x, T y, T z) {
		vec(0, 0) = x;
		vec(1, 0) = y;
		vec(2, 0) = z;
	}
	template <typename T> bool operator()(T const* const* parameter, T* residual) const {
		Eigen::Matrix<T, 3, 1> p_i;
		initailize<T>(p_i, T(px_), T(py_), parameter[0][0]);
		//p_i << T(px_), T(py_), parameter[0];
		Eigen::Matrix<T, 3, 1> p_ip1;
		//p_ip1 << T(px_), T(py_), parameter[1];
		initailize<T>(p_ip1, T(px_), T(py_), parameter[1][0]);


		Eigen::Map<Eigen::Matrix<T, 3, 1>> resMat(residual);

		resMat = (p_ip1 - p_ip1)* weight;

		return true;
	}
	static void configCostFunction(ceres::DynamicAutoDiffCostFunction<L2Regularizer>* pCostFunction) {
		for (size_t j = 0; j < 2; j++)
		{
			pCostFunction->AddParameterBlock(1);
		}

		pCostFunction->SetNumResiduals(3);
	}
	double px_, py_;
	double weight = 0.03;
};

/*
	Takes 4 parameter blocks: vertex i, j position for this and next frame
*/
struct IsometricDeformationFieldLossOnlyZ
{
	IsometricDeformationFieldLossOnlyZ(double px, double py, double qx, double qy):
		px_(px), py_(py), qx_(qx), qy_(qy)
	{}
	/*
		Inputs: vertex_i, vertex_i+1, vertex_j, vertex_j+1, 
	*/

	template <typename T> static void initailize(Eigen::Matrix<T, 3, 1>& vec, T x, T y, T z) {
		vec(0, 0) = x;
		vec(1, 0) = y;
		vec(2, 0) = z;
	}
	template <typename T> bool operator()(T const* const* parameter, T* residual) const {
		Eigen::Matrix<T, 3, 1> p_i;
		initailize<T>(p_i, T(px_), T(py_), parameter[0][0]);
		//p_i << T(px_), T(py_), parameter[0];
		Eigen::Matrix<T, 3, 1> p_ip1;
		//p_ip1 << T(px_), T(py_), parameter[1];
		initailize<T>(p_ip1, T(px_), T(py_), parameter[1][0]);

		Eigen::Matrix<T, 3, 1> q_i;
		//q_i << T(qx_), T(qy_), parameter[2];
		initailize<T>(q_i, T(qx_), T(qy_), parameter[2][0]);

		Eigen::Matrix<T, 3, 1> q_ip1;
		//q_ip1 << T(qx_), T(qy_), parameter[3];
		initailize<T>(q_ip1, T(qx_), T(qy_), parameter[3][0]);

		*residual = ((p_ip1 - p_i - q_ip1 + q_i).transpose() * (q_i - p_i))(0);

		return true;
	}

	static void configCostFunction(ceres::DynamicAutoDiffCostFunction<IsometricDeformationFieldLossOnlyZ>* pCostFunction ) {
		for (size_t j = 0; j < 4; j++)
		{
			pCostFunction->AddParameterBlock(1);
		}

		pCostFunction->SetNumResiduals(1);
	}

	double px_, py_, qx_, qy_;
};

//std::recursive_mutex lock;
void parseArgs(int argc, char** argv, Config& config) {
	std::string modeStr;
	cxxopts::Options options(argv[0], "App generate correspondences and triangles from label files.");

	options.add_options()
		("inExt", "[Optional]Ext name for input files.", cxxopts::value<std::string>(config.inExtName))
		("b,beginFrame", "[Optional]Frame id to start process with.", cxxopts::value<int>(config.startFrameId))
		("e,endFrame", "[Optional]Last frame id to process with.", cxxopts::value<int>(config.endFrameId))
		;

	try
	{
		auto result = options.parse(argc, argv);
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << "\n";
		std::cout << options.help();
		exit(1);
	}

}

int main(int argc, char** argv) {
	clock_t timeStart = clock();
	if (argc < 3)
	{
		std::cout << "Need at least 2 parameters: [inFolder] [outFolder]" << "\n";
		return -1;
	}
	//M mesh;
	////A sample from MeshLab
	//clock_t time1 = clock();
	//mesh.read_ply("F:/WorkingCopy2/2020_11_26_SMPLSHFit/Fit/Lada_Stand/AA00008274.ply");
	//time1 = clock() - time1;
	//printf("Time consumption in loading: %d\n", time1);
	//getchar();
	//

	////Usage of MVIterator, to iterate all vertices in mesh
	//for (M::EPtr pE : It::MEIterator(&mesh))
	//{
	//	printf("pE->id(): %d, edge vertices: %u to %u\n", pE->index(), pE->halfedge()->source()->index(), pE->halfedge()->target()->index());

	//	////Usage of VFIterator, to iterate all faces around vertex
	//	//for (M::FPtr pF : It::MEIterator(pV))
	//	//{
	//	//	printf("\tpF->id: %d : ", pF->id());
	//	//	//Usage of FVIterator, to iterate all vertex around face
	//	//	for (M::VPtr pVF : It::FVIterator(pF))
	//	//	{
	//	//		printf("%d, ", pVF->id());
	//	//	}
	//	//	printf("\n");
	//	//}
	//}

	Config cfg;
	parseArgs(argc, argv, cfg);


	// Step 1: load data: vertices and connections
	// Load vertices
	std::string inFolder = argv[1];
	std::string outFolder = argv[2];

	AC::VecStr inFilesAll;
	AC::IO::getFilesWithExt(inFolder, cfg.inExtName, inFilesAll);
	if (cfg.endFrameId == -1) {
		cfg.endFrameId = inFilesAll.size();
	}
	AC::VecStr inFiles(inFilesAll.begin() + cfg.startFrameId, inFilesAll.begin() + cfg.endFrameId);
	if (inFiles.size() == 0)
	{
		std::cout << "No files found in: " << inFolder << "\n";
		return -1;
	}

	std::vector<M*> meshes;

	for (std::string file : inFiles) {
		M* pMesh = new M;
		if (cfg.inExtName == "ply") {
			pMesh->read_ply(file.c_str());
		}
		else if (cfg.inExtName == "obj") {
			pMesh->read_obj(file.c_str());
		}
		meshes.push_back(pMesh);
	}

	std::vector<std::array<int, 2>> edges;
	for (M::EPtr pE : It::MEIterator(meshes[0]))
	{
		//printf("pE->id(): %d, edge vertices: %u to %u\n", pE->index(), pE->halfedge()->source()->index(), pE->halfedge()->target()->index());
		std::array<int, 2> edge = { pE->halfedge()->source()->index(), pE->halfedge()->target()->index() };
		edges.push_back(edge);
	}

	// make parameter blocks
	std::vector<std::vector<double*>> parameterBlocks;
	for (size_t iFrame = 0; iFrame < meshes.size(); ++iFrame)
	{
		parameterBlocks.push_back(std::vector<double*>());
		for (size_t iV = 0; iV < meshes[iFrame]->numVertices(); ++iV) {
			double* pVBlock = (meshes[iFrame]->vertices())[iV].point().ptr();
			parameterBlocks.back().push_back(pVBlock);
		}
	}

	// Prepare the solver!
	Solver::Options options;
	//options.trust_region_strategy_type = ceres::DOGLEG;
	//options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	//options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
	//options.linear_solver_type = ceres::SPARSE_SCHUR;
	options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
	options.num_threads = 6;
	options.max_num_iterations = cfg.maxIteration;
	options.function_tolerance = cfg.functionTolerance;

	options.minimizer_progress_to_stdout = true;
	Solver::Summary summary;
	Problem problem;

	// the isometric penalty loss
	// the are (nFrames-1)x(nEdge) losses
	for (size_t iFrame = 0; iFrame < meshes.size()-1; ++iFrame)
	{
		for (size_t iE = 0; iE < edges.size(); iE++)
		{
			int iSV = edges[iE][0];
			int iTV = edges[iE][1];

			IsometricDeformationFieldLossOnlyZ* pIsoFunctor = new IsometricDeformationFieldLossOnlyZ(parameterBlocks[iFrame][iSV][0], parameterBlocks[iFrame][iSV][1],
				parameterBlocks[iFrame][iTV][0], parameterBlocks[iFrame][iTV][1]);

			ceres::DynamicAutoDiffCostFunction<IsometricDeformationFieldLossOnlyZ>* cost_function_iso =
				new ceres::DynamicAutoDiffCostFunction<IsometricDeformationFieldLossOnlyZ>(pIsoFunctor);

			std::vector<double*> blocks = { parameterBlocks[iFrame][iSV] + 2, parameterBlocks[iFrame+1][iSV] + 2, parameterBlocks[iFrame][iTV] + 2, parameterBlocks[iFrame+1][iTV] + 2 };

			IsometricDeformationFieldLossOnlyZ::configCostFunction(cost_function_iso);

			problem.AddResidualBlock(cost_function_iso, NULL, blocks);

		}
	}

	// Regularizer on deformation field
	//for (size_t iFrame = 0; iFrame < meshes.size() - 1; ++iFrame)
	//{
	//	for (size_t iV = 0; iV < meshes[iFrame]->numVertices(); ++iV) {
	//		/*double* pVBlock = (meshes[iFrame]->vertices())[iV].point().ptr();
	//		parameterBlocks.back().push_back(pVBlock);*/
	//		L2Regularizer* regularizer = new L2Regularizer(parameterBlocks[iFrame][iV][0], parameterBlocks[iFrame + 1][iV][1]);

	//		ceres::DynamicAutoDiffCostFunction<L2Regularizer>* cost_function_l2 =
	//			new ceres::DynamicAutoDiffCostFunction<L2Regularizer>(regularizer);

	//		L2Regularizer::configCostFunction(cost_function_l2);
	//		std::vector<double*> blocks = { parameterBlocks[iFrame][iV] + 2, parameterBlocks[iFrame + 1][iV] + 2};

	//		problem.AddResidualBlock(cost_function_l2, NULL, blocks);

	//	}
	//}

	// set the first and last frame to constant
	for (size_t iV = 0; iV < meshes[0]->numVertices(); ++iV)
	{
		problem.SetParameterBlockConstant(parameterBlocks[0][iV] + 2);
		problem.SetParameterBlockConstant(parameterBlocks[meshes.size() - 1][iV] + 2);
	}


	Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	for (size_t iFrame = 0; iFrame < meshes.size(); ++iFrame)
	{
		AC::IO::FileParts fp(inFiles[iFrame]);
		meshes[iFrame]->write_ply((outFolder + "/" + fp.name + fp.ext).c_str());
	}

}