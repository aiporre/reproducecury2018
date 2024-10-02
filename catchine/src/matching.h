
template <typename TYPE>
int matching(const char *DataFile,const char *ResFile, double RegWeight, int verbosemode) {

	double tic = (double) clock() / (double) CLOCKS_PER_SEC;

	// set the matching class
	Match < TYPE, __Dim__ > Test;
	ifstream f;
	f.open(DataFile);
	try 
	{
		Test.Init(f, RegWeight, verbosemode);
	} 
	catch (int ret) 
	{
		f.close();
		return ret;
	}
	f.close();

	// Optimization
	Test.DoOptimize();

	double toc = (double) clock() / (double) CLOCKS_PER_SEC;
	if (verbosemode) cout << toc - tic << " seconds" << endl;

	// output
	if (verbosemode) cout << "Writing results to file " << ResFile << endl;
	ofstream fout;
	fout.open(ResFile);
	Test.Write(fout);
	fout.close();

	return 0;
}
