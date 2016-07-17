#include "DenoisingFacade.h"

#include "Algorithms/BilateralMeshDenoising.h"
#include "Algorithms/NonIterativeFeaturePreservingMeshFiltering.h"
#include "Algorithms/FastAndEffectiveFeaturePreservingMeshDenoising.h"
#include "Algorithms/BilateralNormalFilteringForMeshDenoising.h"
#include "Algorithms/MeshDenoisingViaL0Minimization.h"
#include "Algorithms/GuidedMeshNormalFiltering.h"

using std::pair;
//enum AlgorithmsType {
//	kNon, kNoise,
//	kBilateralMeshDenoising, kNonIterativeFeaturePreservingMeshFiltering,
//	kFastAndEffectiveFeaturePreservingMeshDenoising, kBilateralNormalFilteringForMeshDenoising,
//	kMeshDenoisingViaL0Minimization, kGuidedMeshNormalFiltering
DenoisingFacade::DenoisingFacade()
{
	noise_ = NULL;
	mesh_denoise_base_ = NULL;
	algorithms_type_ = AlgorithmsType::kNon;

	algorithms_dictionary_.insert(std::make_pair("Noise", AlgorithmsType::kNoise));
	algorithms_dictionary_.insert(std::make_pair("BilateralMeshDenoising", AlgorithmsType::kBilateralMeshDenoising));
	algorithms_dictionary_.insert(std::make_pair("NonIterativeFeaturePreservingMeshFiltering", AlgorithmsType::kNonIterativeFeaturePreservingMeshFiltering));
	algorithms_dictionary_.insert(std::make_pair("FastAndEffectiveFeaturePreservingMeshDenoising", AlgorithmsType::kFastAndEffectiveFeaturePreservingMeshDenoising));
	algorithms_dictionary_.insert(std::make_pair("BilateralNormalFilteringForMeshDenoising", AlgorithmsType::kBilateralNormalFilteringForMeshDenoising));
	algorithms_dictionary_.insert(std::make_pair("MeshDenoisingViaL0Minimization", AlgorithmsType::kMeshDenoisingViaL0Minimization));
	algorithms_dictionary_.insert(std::make_pair("GuidedMeshNormalFiltering", AlgorithmsType::kGuidedMeshNormalFiltering));
	

	noise_type_dictionary_.insert(std::make_pair("Gaussian", Noise::NoiseType::kGaussian));
	noise_type_dictionary_.insert(std::make_pair("Impulsive", Noise::NoiseType::kImpulsive));
}


DenoisingFacade::~DenoisingFacade()
{
	if (noise_ != NULL)
		delete noise_;

	if (mesh_denoise_base_ != NULL)
		delete mesh_denoise_base_;
}


void DenoisingFacade::initAlgorithm(DataManager *_data_manager, ParameterSet *_parameter_set)
{
	switch (algorithms_type_) {
	case kNoise:
		noise_ = new Noise(_data_manager, _parameter_set);
		break;
	case kBilateralMeshDenoising:
		mesh_denoise_base_ = new BilateralMeshDenoising(_data_manager, _parameter_set);
		break;
	case kNonIterativeFeaturePreservingMeshFiltering:
		mesh_denoise_base_ = new NonIterativeFeaturePreservingMeshFiltering(_data_manager, _parameter_set);
		break;
	case kFastAndEffectiveFeaturePreservingMeshDenoising:
		mesh_denoise_base_ = new FastAndEffectiveFeaturePreservingMeshDenoising(_data_manager, _parameter_set);
		break;
	case kBilateralNormalFilteringForMeshDenoising:
		mesh_denoise_base_ = new BilateralNormalFilteringForMeshDenoising(_data_manager, _parameter_set);
		break;
	case kMeshDenoisingViaL0Minimization:
		mesh_denoise_base_ = new MeshDenoisingViaL0Minimization(_data_manager, _parameter_set);
		break;
	case kGuidedMeshNormalFiltering:
		mesh_denoise_base_ = new GuidedMeshNormalFiltering(_data_manager, _parameter_set);
		break;
	default:
		break;
	}
}

void DenoisingFacade::run()
{
	if (algorithms_type_ == kNoise)
		noise_->addNoise();
	else
		mesh_denoise_base_->denoise();
}
void DenoisingFacade::setAlgorithmType(const string& type)
{
	algorithms_type_ = algorithms_dictionary_.at(type);
}
Noise::NoiseType DenoisingFacade::getNoiseType(const string& type)
{
	return noise_type_dictionary_.at(type);
}