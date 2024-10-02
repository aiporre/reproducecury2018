
template <typename TYPE>
int matching(const char* DataFile,const char* ResFile, int niters, double StepSize, double RegWeight, int useoptim, double breakratio, int loopbreak, int verbosemode)
{
	
    double tic = (double) clock () / (double) CLOCKS_PER_SEC;
	
    // set the matching class
	Match < TYPE, __Dim__ > Test;
	ifstream f;
    f.open(DataFile);
	try
	{
		Test.Init(f, RegWeight, verbosemode);
	}
	catch(int ret)
	{
		f.close();
		return ret;
    }
	f.close();
	
    // Optimization
    switch(useoptim)
	{
		case 0:
			Test.DoFixedGradDescent(StepSize,niters);
			break;
        case 1:
			Test.DoAdaptGradDescent(StepSize, niters, breakratio, loopbreak);
			break;	
        case 2:
			Test.DoLBFGS();
			break;
#ifdef DLIB_CATCHINE
		case 3:
			Test.DoLBFGS_dlib();
			break;
		case 4:
			Test.DoBFGS_dlib();
			break;
		case 5:
			Test.DoCG_dlib();
			break;
		case 6:
			Test.DoLBFGS_approxder_dlib();
			break;
#endif
    }
	
	
    double toc = (double) clock () / (double) CLOCKS_PER_SEC;
    if(verbosemode) cout << toc - tic << " seconds" << endl;
	
    // output
    if(verbosemode) cout << "Writing results to file " << ResFile << endl;
    ofstream fout;
    fout.open (ResFile);
    Test.Write(fout);
    fout.close ();
	
    return 0;
}
