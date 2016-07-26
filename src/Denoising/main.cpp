
#include <QtCore/QCoreApplication>
#include <QCommandLineParser>
#include <QFileInfo> 
#include <iostream>
#include <QElapsedTimer>
#include "DenoisingFacade.h"

using std::cout; using std::endl;

int main(int argc, char *argv[])
{
	QCoreApplication a(argc, argv);	
	QCoreApplication::setApplicationVersion("0.1 2016-07-13");

	QCommandLineParser parser;
	parser.setApplicationDescription("Is this the help?");
	parser.addHelpOption();//Adds the help option(-h, --help and -? on Windows)
	parser.addVersionOption();//-v / --version option, which displays the version string of the application.
	parser.addPositionalArgument("source", QCoreApplication::translate("main", "Source file to denoise."));
	parser.addPositionalArgument("destination", QCoreApplication::translate("main", "Destination directory."));
	parser.addPositionalArgument("algorithmType", QCoreApplication::translate("main", "Algorithm type: \
																					  Noise, \
																					  BilateralMeshDenoising, \
																					  NonIterativeFeaturePreservingMeshFiltering, \
																					  FastAndEffectiveFeaturePreservingMeshDenoising, \
																					  BilateralNormalFilteringForMeshDenoising, \
																					  MeshDenoisingViaL0Minimization, \
																					  GuidedMeshNormalFiltering, \
																					  ShortestPropagationMeshFiltering."));
	
	// options for noise
	// todo relation between noise level and impulsive level
	// add more optioins: Impulsive level, Noise direction
	QCommandLineOption noiseTypeOption(QStringList() << "noiseType", 
		"Noise type: Gaussian, Impulsive.", "Noise type", "Gaussian");
	parser.addOption(noiseTypeOption);
	QCommandLineOption noiseLevelOption(QStringList() << "noiseLevel", 
		"Noise Level: between 0 and 1.", "Noise level", "0.5");
	parser.addOption(noiseLevelOption);
	

	// options for GuidedMeshNormalFiltering && ShortestPropagationMeshFiltering
	QCommandLineOption FaceNeighborTypeOption(QStringList() << "faceNeighborType",
		"The type of the neighbor of the face: geometrical (radius based: 2) or topological (vertex: 0, edge: 1 or face-ring based: 3)", "Face neighbor type", "geometrical");
	parser.addOption(FaceNeighborTypeOption);
	QCommandLineOption FaceDistOption(QStringList() << "FaceDist",
		"Multiple(* avg face dis.) => Radius for search geometrical neighbor of the face; or it means topological distance", "Face distant");
	parser.addOption(FaceDistOption);

	QCommandLineOption SigmaSOption(QStringList() << "SigmaS", "sigma_s", "Multiple(* sigma_s)");
	parser.addOption(SigmaSOption);
	QCommandLineOption SigmaROption(QStringList() << "SigmaR", "sigma_r", "sigma_r");
	parser.addOption(SigmaROption);

	QCommandLineOption NormalIterNumOption(QStringList() << "NormalIterNum", "Normal Iteration Number", "(Local)Normal Iteration Num.");
	parser.addOption(NormalIterNumOption);
	QCommandLineOption VertexIterNumOption(QStringList() << "VertexIterNum", "Vertex Iteration Number", "Vertex Iteration Num.");
	parser.addOption(VertexIterNumOption);

	// Process(处理) the actual command line arguments given by the user
	parser.process(a);


	////////////////////////////////////////////////////////////
	////  load mesh         ////////////////////////////////////
	////////////////////////////////////////////////////////////
	const QStringList args = parser.positionalArguments();
	// source is args.at(0), destination is args.at(1)
	QFileInfo finfo(args.at(0));
	QString inFilePath(finfo.filePath()); // ("..\\..\\models\\Fandisk0.3\\Original.obj");
	QString outDir(args.at(1));

	DataManager dm;
	if (!dm.ImportMeshFromFile( inFilePath.toStdString() ) )//转化成标准的字符串
	{
		cout << "Loading mesh " << inFilePath.toStdString() << " failed." << endl;;
		return -1;
	}
	else
	{
		cout << "Loading mesh " << inFilePath.toStdString() << ", num of faces: " << dm.getOriginalMesh().n_faces() << endl;
	}
	

	////////////////////////////////////////////////////////////
	////  facade   （外观外形 ////////////////////////////////////
	////////////////////////////////////////////////////////////

	DenoisingFacade df;
	df.setAlgorithmType( args.at(2).toStdString() );
	ParameterSet params;
	df.initAlgorithm(&dm, &params);//默认参数设定

	////////////////////////////////////////////////////////////
	////  load parameters         ////////////////////////////////////
	////////////////////////////////////////////////////////////

	switch (df.getAlgorithmType())
	{
	case DenoisingFacade::kNoise:
	{
		//更改参数
		QString noiseType = parser.value(noiseTypeOption);
		if (!noiseType.isEmpty())
			params.setStringListIndex("Noise type", df.getNoiseType(noiseType.toStdString()));
		QString noiselevel = parser.value(noiseLevelOption);
		if (!noiselevel.isEmpty())
		{
			if (Noise::NoiseType::kGaussian == df.getNoiseType(noiseType.toStdString()))
			{
				params.setValue("Noise level", noiselevel.toDouble());
				params.setValue("Impulsive level", 0.0);
			}
			else
			{
				params.setValue("Noise level", 0.5);//默认的噪声水平
				params.setValue("Impulsive level", noiselevel.toDouble());
			}
		}
	}
		break;
	default:
		QString FaceDist = parser.value(FaceDistOption);
		if (!FaceDist.isEmpty())
			params.setValue("Multiple(* avg face dis.)", FaceDist.toDouble());

		QString SigmaS = parser.value(SigmaSOption);
		if (!SigmaS.isEmpty())
			params.setValue("Multiple(* sigma_s)", SigmaS.toDouble());

		QString SigmaR = parser.value(SigmaROption);
		if (!SigmaR.isEmpty())
			params.setValue("sigma_r", SigmaR.toDouble());

		QString NormalIterNum = parser.value(NormalIterNumOption);
		if (!NormalIterNum.isEmpty())
			params.setValue("(Local)Normal Iteration Num.", NormalIterNum.toInt());

		QString VertexIterNum = parser.value(VertexIterNumOption);
		if (!VertexIterNum.isEmpty())
			params.setValue("Vertex Iteration Num.", VertexIterNum.toInt());
	}


	////////////////////////////////////////////////////////////
	////       compute  and  output    //////////////////////////////////
	////////////////////////////////////////////////////////////	
	QElapsedTimer timer;
	timer.start();
	df.run();
	std::cout << "time elapsed: " << timer.elapsed() << std::endl;

	if (DenoisingFacade::kNoise == df.getAlgorithmType())
		dm.MeshToNoisyMesh();
	else
		dm.MeshToDenoisedMesh();	

	QString outFilePath(outDir + '/' + finfo.baseName() + "_result." + finfo.suffix());
	if (!dm.ExportMeshToFile(outFilePath.toStdString()) )
	{
		cout << "Writing mesh " << outFilePath.toStdString() << " failed." << endl;;
		return -1;
	}
	else
		cout << "Writing mesh " << outFilePath.toStdString() << " successful." << endl;

	return 0;
	//return a.exec();
}
