
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
	parser.addPositionalArgument("algorithmType", QCoreApplication::translate("main", "Algorithm type: Noise, BilateralMeshDenoising, NonIterativeFeaturePreservingMeshFiltering, FastAndEffectiveFeaturePreservingMeshDenoising, BilateralNormalFilteringForMeshDenoising, MeshDenoisingViaL0Minimization, GuidedMeshNormalFiltering."));
	

	QCommandLineOption noiseTypeOption(QStringList() << "noiseType", QCoreApplication::translate("main", "Noise type: Gaussian, Impulsive."), "Noise type" );
	parser.addOption(noiseTypeOption);
	QCommandLineOption noiseLevelOption(QStringList() << "noiseLevel", QCoreApplication::translate("main", "Noise Level: between 0 and 1."), "Noise level");
	parser.addOption(noiseLevelOption);


	QCommandLineOption GuidedFilterFaceDistOption(QStringList() << "GuidedFilterFaceDist", QCoreApplication::translate("main", "how much around a face."), "Multiple(* avg face dis.)");
	parser.addOption(GuidedFilterFaceDistOption);
	QCommandLineOption GuidedFilterSigmaSOption(QStringList() << "GuidedFilterSigmaS", QCoreApplication::translate("main", "sigma_s"), "Multiple(* sigma_s)");
	parser.addOption(GuidedFilterSigmaSOption);
	QCommandLineOption GuidedFilterSigmaROption(QStringList() << "GuidedFilterSigmaR", QCoreApplication::translate("main", "sigma_r"), "sigma_r");
	parser.addOption(GuidedFilterSigmaROption);
	QCommandLineOption GuidedFilterNormalIterNumOption(QStringList() << "GuidedFilterNormalIterNum", QCoreApplication::translate("main", "Normal Iteration Number"), "(Local)Normal Iteration Num.");
	parser.addOption(GuidedFilterNormalIterNumOption);
	QCommandLineOption GuidedFilterVertexIterNumOption(QStringList() << "GuidedFilterVertexIterNum", QCoreApplication::translate("main", "Vertex Iteration Number"), "Vertex Iteration Num.");
	parser.addOption(GuidedFilterVertexIterNumOption);

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
	df.initAlgorithm(&dm, &params);//默认参数设定

	////////////////////////////////////////////////////////////
	////  load parameters         ////////////////////////////////////
	////////////////////////////////////////////////////////////

	if (args.at(2).toStdString() == "Noise")
	{
		//更改参数
		QString noiseType = parser.value(noiseTypeOption);
		if (!noiseType.isEmpty())
			params.setStringListIndex("Noise type", df.getNoiseType(noiseType.toStdString()));
		QString noiselevel = parser.value(noiseLevelOption);
		if (!noiselevel.isEmpty())
		{
			if (noiseType.toStdString() == "Gaussian")
			{
				params.setValue("Noise level", df.getParaValue(noiselevel.toStdString()));
				params.setValue("Impulsive level", 0.0);
			}
				
			if (noiseType.toStdString() == "Impulsive")
			{	
				params.setValue("Noise level", 0.5);//默认的噪声水平
				params.setValue("Impulsive level", df.getParaValue(noiselevel.toStdString()));
				
			}
				
		}
	}
	if (args.at(2).toStdString() == "GuidedMeshNormalFiltering")
	{
		QString FaceDist = parser.value(GuidedFilterFaceDistOption);
		if (!FaceDist.isEmpty())
			params.setValue("Multiple(* avg face dis.)",df.getParaValue(FaceDist.toStdString()));
		QString SigmaS = parser.value(GuidedFilterSigmaSOption);
		if (!SigmaS.isEmpty())
			params.setValue("Multiple(* sigma_s)", df.getParaValue(SigmaS.toStdString()));
		QString SigmaR = parser.value(GuidedFilterSigmaROption);
		if (!SigmaR.isEmpty())
			params.setValue("sigma_r", df.getParaValue(SigmaR.toStdString()));
		QString NormalIterNum = parser.value(GuidedFilterNormalIterNumOption);
		if (!NormalIterNum.isEmpty())
			params.setValue("(Local)Normal Iteration Num", (int)df.getParaValue(NormalIterNum.toStdString()));
		QString VertexIterNum = parser.value(GuidedFilterVertexIterNumOption);
		if (!VertexIterNum.isEmpty())
			params.setValue("Vertex Iteration Num", (int)df.getParaValue(VertexIterNum.toStdString()));
	}


		////////////////////////////////////////////////////////////
		////       compute  and  output    //////////////////////////////////
		////////////////////////////////////////////////////////////
		
	df.run();
	if (args.at(2).toStdString() == "Noise")
	{
		dm.MeshToNoisyMesh();
	}
	else
	{
		dm.MeshToDenoisedMesh();
	}

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
