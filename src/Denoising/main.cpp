
#include <QtCore/QCoreApplication>
#include <QCommandLineParser>
#include <QFileInfo> 
#include <iostream>
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

	parser.addPositionalArgument("algorithmType", QCoreApplication::translate("main", "Algorithm type: Noise, ..."));
	

	QCommandLineOption noiseTypeOption(QStringList() << "noiseType", QCoreApplication::translate("main", "Noise type: Gaussian, Impulsive."), "Noise type" );
	parser.addOption(noiseTypeOption);


	QCommandLineOption AlgorithmsTypeOption(QStringList() << "AlgorithmsType", 
		QCoreApplication::translate("main", "Algorithms type: BilateralMeshDenoising, NonIterativeFeaturePreservingMeshFiltering, FastAndEffectiveFeaturePreservingMeshDenoising, BilateralNormalFilteringForMeshDenoising, MeshDenoisingViaL0Minimization, GuidedMeshNormalFiltering."),
		"Algorithms type");
	parser.addOption(AlgorithmsTypeOption);

	// Process(处理) the actual command line arguments given by the user
	parser.process(a);


	const QStringList args = parser.positionalArguments();

	// source is args.at(0), destination is args.at(1)
	QFileInfo finfo( args.at(0) );	

	QString inFilePath( finfo.filePath() ); // ("..\\..\\models\\Fandisk0.3\\Original.obj");

	QString outDir( args.at(1) );

	
	
	////////////////////////////////////////////////////////////
	////  load mesh         ////////////////////////////////////
	////////////////////////////////////////////////////////////
	DataManager dm;
	if (!dm.ImportMeshFromFile( inFilePath.toStdString() ) )//转化成标准的字符串
	{
		cout << "Loading mesh " << inFilePath.toStdString() << " failed." << endl;;
		return -1;
	}
	else
		cout << "Loading mesh " << inFilePath.toStdString() << " successful." << endl;
	

	////////////////////////////////////////////////////////////
	////  facade   （外观外形 ////////////////////////////////////
	////////////////////////////////////////////////////////////
	DenoisingFacade df;
	df.setAlgorithmType( args.at(2).toStdString() );
	ParameterSet params;
	df.initAlgorithm(&dm, &params);//默认初始化有些参数设定

	////////////////////////////////////////////////////////////
	////  load parameters         ////////////////////////////////////
	////////////////////////////////////////////////////////////
	QString noiseType = parser.value(noiseTypeOption);
	if (!noiseType.isEmpty())
		params.setStringListIndex("Noise type", df.getNoiseType(noiseType.toStdString()));//参数不够可以增加函数

	QString AlgorithmsType = parser.value(AlgorithmsTypeOption);
	if (!noiseType.isEmpty())
		params.setStringListIndex("Algorithms type", df.getNoiseType(AlgorithmsType.toStdString()));//参数不够可以增加函数

		
	////////////////////////////////////////////////////////////
	////  compute and output //////////////////////////////////
	////////////////////////////////////////////////////////////
	df.run();
	if (!noiseType.isEmpty())
		dm.MeshToNoisyMesh();
	if (!noiseType.isEmpty())
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
