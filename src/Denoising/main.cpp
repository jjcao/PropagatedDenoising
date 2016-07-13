
#include <QtCore/QCoreApplication>
#include <QCommandLineParser>
#include <QFileInfo> 
#include <iostream>
#include "DenoisingFacade.h"

using std::cout; using std::endl;

int main(int argc, char *argv[])
{
	QCoreApplication a(argc, argv);	
	QCommandLineParser parser;
	QCoreApplication::setApplicationVersion("0.1 2016-07-13");
	parser.setApplicationDescription("Is this the help?");
	parser.addHelpOption();//Adds the help option(-h, --help and -? on Windows)
	parser.addVersionOption();//-v / --version option, which displays the version string of the application.
	parser.addPositionalArgument("source", QCoreApplication::translate("main", "Source file to denoise."));
	parser.addPositionalArgument("destination", QCoreApplication::translate("main", "Destination directory."));
	parser.addPositionalArgument("algorithmType", QCoreApplication::translate("main", "Algorithm type: Noise, ..."));
	

	QCommandLineOption noiseTypeOption(QStringList() << "noiseType",
		QCoreApplication::translate("main", "Noise type: Gaussian, Impulsive."), "Noise type" );
	parser.addOption(noiseTypeOption);


	parser.process(a);
	const QStringList args = parser.positionalArguments();

	QFileInfo finfo( args.at(0) );	
	QString inFilePath( finfo.filePath() ); // ("..\\..\\models\\Fandisk0.3\\Original.obj");
	QString outDir( args.at(1) );

	
	
	////////////////////////////////////////////////////////////
	////  load mesh         ////////////////////////////////////
	////////////////////////////////////////////////////////////
	DataManager dm;
	if (!dm.ImportMeshFromFile( inFilePath.toStdString() ) )
	{
		cout << "Loading mesh " << inFilePath.toStdString() << " failed." << endl;;
		return -1;
	}
	else
		cout << "Loading mesh " << inFilePath.toStdString() << " successful." << endl;
	

	////////////////////////////////////////////////////////////
	////  facade         ////////////////////////////////////
	////////////////////////////////////////////////////////////
	DenoisingFacade df;
	df.setAlgorithmType( args.at(2).toStdString() );
	ParameterSet params;
	df.initAlgorithm(&dm, &params);

	////////////////////////////////////////////////////////////
	////  load parameters         ////////////////////////////////////
	////////////////////////////////////////////////////////////
	QString noiseType = parser.value(noiseTypeOption);
	if (!noiseType.isEmpty())
		params.setStringListIndex("Noise type", df.getNoiseType(noiseType.toStdString()));

		
	////////////////////////////////////////////////////////////
	////  compute and output //////////////////////////////////
	////////////////////////////////////////////////////////////
	df.run();
	dm.MeshToNoisyMesh();
	QString outFilePath(outDir + '/' + finfo.baseName() + '.' + finfo.suffix());
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
