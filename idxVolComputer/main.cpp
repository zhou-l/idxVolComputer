#include "IdxVolComputer.h"
#include "global.h"
#include <QtWidgets/QApplication>


int main(int argc, char* argv[])
{
	QApplication::setAttribute(Qt::AA_ShareOpenGLContexts);
	QApplication a(argc, argv);

	// For now, set file names in the argument list
	//if (argc < 3 || argc > 4)
	//{
	//	fprintf(stderr, "\vol_cont_scatterplot [NRRD_input_file_name] [scatterplot_output_filename]\nOr: \n");
	//	fprintf(stderr, "\vol_cont_scatterplot [NRRD_input_file_name] [NRRD_input_file_name2] [scatterplot_output_filename]\n\n");
	//	return -1;
	//}

	// global variable initialization
	global_init();
	IdxVolComputer w;

	QStringList  dataList;

	for (int c = 1; c < argc; c++)
		dataList.push_back(argv[c]);
	if (!dataList.empty()) {
		if (!w.loadFiles(dataList))
		{
			w.cleanData();
			return -1;
		}
		else
			w.drawPlot();
	}



	w.show();

	return a.exec();
}
