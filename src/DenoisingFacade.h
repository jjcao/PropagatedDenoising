#pragma once

#include <string>
#include <map>
#include "Algorithms/Noise.h"
#include "Algorithms/MeshDenoisingBase.h"

using std::string;
using std::map;

class DenoisingFacade
{
public:
	DenoisingFacade();
	~DenoisingFacade();

	enum AlgorithmsType {
		kNon,                                            kNoise, 
		kBilateralMeshDenoising,                         kNonIterativeFeaturePreservingMeshFiltering,
		kFastAndEffectiveFeaturePreservingMeshDenoising, kBilateralNormalFilteringForMeshDenoising,
		kMeshDenoisingViaL0Minimization,                 kGuidedMeshNormalFiltering
	};

	void setAlgorithmType(const string& type);
	void initAlgorithm(DataManager *_data_manager, ParameterSet *_parameter_set);
	Noise::NoiseType getNoiseType(const string& type);
	void run();

private:
	AlgorithmsType algorithms_type_;
	map<string, AlgorithmsType> algorithms_dictionary_;
	map<string, Noise::NoiseType> noise_type_dictionary_;
	Noise *noise_;
	MeshDenoisingBase *mesh_denoise_base_;
};

