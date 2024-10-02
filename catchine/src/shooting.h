
template <typename TYPE>
int shooting(char* DataFile, char* ResFile)
{
    LargeDef_InitParam<TYPE,__Dim__> *Ev = NULL;

    ifstream f;
    f.open(DataFile);
    string Tag;
    f >> Tag;
    if(!Tag.compare("Evol"))
    {
        f >> Tag;
        if(!Tag.compare("LargeDef_InitParam"))
            Ev = new LargeDef_InitParam < TYPE, __Dim__ > (f);
		else
        {
            cout << "Should be LargeDef_InitParam; exiting" << endl;
            f.close();
			return -1;
        }
    }
    else
    {
        cout << "error : first line of "<< ResFile << " should read : Evol" <<endl;
        f.close();
		return -1;
    }
    f.close();


    // shooting
    cout << "Starting shooting.." << endl;
    double tic = (double) clock () / (double) CLOCKS_PER_SEC;
    Ev->AllShooting();
    double toc = (double) clock () / (double) CLOCKS_PER_SEC;
    cout << toc - tic << " seconds for Shooting" << endl;

    // output
    cout << "Writing results to file " << ResFile << endl;
    ofstream fout;
    fout.open (ResFile);
    fout << "Evol" << endl;
    Ev->Write(fout);
    fout.close ();
	delete(Ev);
    return 0;
}
