
template <typename TYPE>
int flow(const char* ResFile, const char* PointsFileIn, const char* PointsFileOut, int tdeb, int tfin)
{
    typedef TinyVector<TYPE,__Dim__> Vect;
    typedef Array<Vect,1> ArrVect;
    typedef Array<ArrVect,1> ArrArrVect;
    // set the Deformation object
    Evol < TYPE, __Dim__ > *Ev = NULL;
    ifstream f;
    f.open(ResFile);
    string Tag;
    f >> Tag;
    if(!Tag.compare("Evol"))
    {
        f >> Tag;
        if(!Tag.compare("LargeDef_InitParam"))
            Ev = new LargeDef_InitParam < TYPE, __Dim__ > (f);
        else if(!Tag.compare("LargeDef"))
            Ev = new LargeDef < TYPE, __Dim__ > (f);
    }
    else
    {
        cout << "error : first line of "<< ResFile << " should read : Evol" <<endl;
        f.close();
		return -1;
    }
    f.close();

    // reading the Points to flow
    ArrVect P;
    f.open(PointsFileIn);
    ReadArr(P,f);
    f.close();

    ArrVect PhiP(Range(1,P.rows()));
    double tic, toc;
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    Ev->Flow(P,PhiP,tdeb,tfin);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
    cout << toc - tic << " seconds for Flow" << endl;

    // output
    ofstream fout;
    fout.open (PointsFileOut);
    WriteArr(PhiP,fout);
    fout.close ();
	delete(Ev);
    return 0;
}
